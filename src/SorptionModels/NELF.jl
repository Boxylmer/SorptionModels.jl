struct NELF end


"""
    Requires that in the kijmatrix, the polymer is the first index. The remaining indexes must match `penetrants`.
"""

struct NELFModel{BMT, POLYMT, PDT, KSWT} <: SorptionModel
    bulk_model::BMT                 # EOS
    polymer_model::POLYMT           # EOS
    polymer_dry_density::PDT        # number
    ksw_values::KSWT                # vector of values
end
# todo, add =[1] to bulk_phase_mole_fractions as a default constructor where applicable
function bulk_phase_activity(model::NELFModel, temperature, pressure, bulk_phase_mole_fractions)
    μ = activity(
        model.bulk_model, 
        pressure, 
        temperature, 
        bulk_phase_mole_fractions)
    return μ
end

function predict_concentration(model::NELFModel, temperature::Number, pressure::Number, bulk_phase_mole_fractions; units=:cc)
    pressure = pressure < eps() ? eps() : pressure    
    bulk_phase_mole_fractions = ((i) -> (i < eps() ? eps() : i)).(bulk_phase_mole_fractions)  # Courtesy of Clementine (Julia Discord)

    target_activities = bulk_phase_activity(model, temperature, pressure, bulk_phase_mole_fractions)
    penetrant_mass_fraction_initial_guesses = ones(length(bulk_phase_mole_fractions)) * 1e-3
    
    function activity_error(penetrant_mass_fractions)
        polymer_mass_fraction = 1 - sum(penetrant_mass_fractions)
        if polymer_mass_fraction <= 0 || polymer_mass_fraction >= 1
            @warn "Polymer mass fraction was not valid, returning very large error value."
            @show polymer_mass_fraction
            return 1e100
        end
        polymer_phase_mass_fractions = vcat(polymer_mass_fraction, penetrant_mass_fractions)
        polymer_phase_density_after_swelling = calculate_polymer_phase_density(model, pressure, bulk_phase_mole_fractions, polymer_phase_mass_fractions)
        polymer_phase_density_upper_bound = density_upper_bound(model.polymer_model, polymer_phase_mass_fractions)
        if polymer_phase_density_after_swelling > polymer_phase_density_upper_bound 
            return 1e100
        end

        polymer_phase_activities = ρTω_activity(
            model.polymer_model, 
            polymer_phase_density_after_swelling, 
            temperature, 
            polymer_phase_mass_fractions)
        residual_squared = rss(target_activities, polymer_phase_activities[2:end]) 
        return residual_squared
    end
    lower = ones(length(penetrant_mass_fraction_initial_guesses))*eps()
    upper = ones(length(penetrant_mass_fraction_initial_guesses)) .- eps()
    res = Optim.optimize(activity_error, lower, upper, penetrant_mass_fraction_initial_guesses, Fminbox(LBFGS()), Optim.Options(
        allow_f_increases = false,
        x_tol = 1e-5,
        g_tol = 1e-7
    ); autodiff=:forward)

    penetrant_mass_fractions = Optim.minimizer(res)
    polymer_phase_mass_fractions = vcat(1 - sum(penetrant_mass_fractions), penetrant_mass_fractions)

    if units==:frac
        return polymer_phase_mass_fractions
    elseif units==:g
        concs_g_g = polymer_phase_mass_fractions_to_gpen_per_gpol(polymer_phase_mass_fractions)
        return concs_g_g
    elseif units==:cc
        concs_cc_cc = polymer_phase_mass_fractions_to_ccpen_per_ccpol(polymer_phase_mass_fractions, model.polymer_dry_density, molecular_weight(model.bulk_model))
        return concs_cc_cc
    end
end
predict_concentration(model::NELFModel, temperature::Number, pressure::Number; kwargs...) = predict_concentration(model, temperature, pressure, [1]; kwargs...)

function calculate_swelled_polymer_density(model::NELFModel, penetrant_partial_pressures)
    return model.polymer_dry_density * (1 - sum(model.ksw_values .* penetrant_partial_pressures))
end

function calculate_swelled_polymer_density(model::NELFModel, pressure, penetrant_mole_fractions)
    partial_pressures = penetrant_mole_fractions .* pressure
    return calculate_swelled_polymer_density(model, partial_pressures)
end

function calculate_polymer_phase_density(swelled_density::Number, polymer_phase_mass_fractions::AbstractVector{<:Number})
    return swelled_density / polymer_phase_mass_fractions[1]
end

function calculate_polymer_phase_density(model::NELFModel, pressure::Number, bulk_phase_mole_fractions::AbstractVector{<:Number}, polymer_phase_mass_fractions::AbstractVector{<:Number})
    polymer_density_after_swelling = calculate_swelled_polymer_density(model, pressure, bulk_phase_mole_fractions)
    return calculate_polymer_phase_density(polymer_density_after_swelling, polymer_phase_mass_fractions)
end


# functions for fit data to NELF parameters
"""
    fit_model(::NELF, model_choice, isotherms, bulk_phase_characteristic_params, [polymer_molecular_weight])
Find the EOS parameters of a polymer from a vector of `IsothermData`s using the NELF model. 
# Arguments
- `model_choice`: MembraneEOS model to use
- `isotherms`: `Vector` of `IsothermData` structs using the polymer in question. Temperature, pressure, density, and concentration must be provided.
- `bulk_phase_characteristic_params`: Vector of pure characteristic parameter vectors following the same order as the isotherms. 
  - E.g., for Sanchez Lacombe and two input isotherms, `bulk_phase_characteristic_params = [[p★_1, t★_1, ρ★_1, mw_1], [p★_2, t★_2, ρ★_2, mw_2]]`
"""
function fit_model(::NELF, isotherms::AbstractVector{<:IsothermData}, bulk_phase_characteristic_params; 
    polymer_molecular_weight=10000, eos=MembraneEOS.SanchezLacombe())
    
    # this function uses SL, which needs 4 params per component, one of which is already specified (MW)
    
    dualmode_models = [fit_model(DualMode(), isotherm) for isotherm in isotherms]
    
    bulk_phase_models = [SL(params...) for params in bulk_phase_characteristic_params]
    default_ksw_vec = [0]
    densities = polymer_density.(isotherms) # get each isotherm's density in case the user accounted for polymers from different batches
    temperatures = temperature.(isotherms)
    infinite_dilution_pressure = 1e-4 # ???

    function error_function(char_param_vec)
        given_sol = zeros(length(isotherms))
        pred_sol = zeros(length(isotherms))

        for i in eachindex(isotherms)
            char_pressures = [char_param_vec[1], bulk_phase_characteristic_params[i][1]]
            char_temperatures = [char_param_vec[2], bulk_phase_characteristic_params[i][2]]
            char_densities = [char_param_vec[3], bulk_phase_characteristic_params[i][3]]
            molecular_weights = [polymer_molecular_weight, bulk_phase_characteristic_params[i][4]]

            polymer_phase_model = SL(char_pressures, char_temperatures, char_densities, molecular_weights)
            nelf_model = NELFModel(bulk_phase_models[i], polymer_phase_model, densities[i], default_ksw_vec)
            pred_sol[i] = predict_concentration(nelf_model, temperatures[i], infinite_dilution_pressure, [1])[1]
            given_sol[i] = predict_concentration(dualmode_models[i]::DualModeModel, infinite_dilution_pressure::Number)
        end
        err = rss(given_sol, pred_sol)
        return err
    end
    density_lower_bound = maximum(densities)
    lower = [0., 0., density_lower_bound]
    upper = [3000, 3000, 5.]
    res = Optim.optimize(
        error_function, lower, upper, 
        [500, 500, density_lower_bound * 1.2], 
        Fminbox(LBFGS(; m=60, linesearch = Optim.LineSearches.BackTracking())), 
        Optim.Options(; allow_f_increases = false))
    # res = Optim.optimize(error_function, lower, upper, [500, 500, density_lower_bound * 1.2], SAMIN(; rt = 0.1), Optim.Options(iterations=10^6))

    return [Optim.minimizer(res)..., polymer_molecular_weight]
    # work in progress

end

"""
    fit_kij(NELF(), isotherms, bulk_parameters, polymer_parameters; [interpolation_model]=DualMode(), [kij_fit_p_mpa]=1e-4)
Find the best kij and ksw parameters for the NELF model according to the Sanchez Lacombe EOS. 
# Arguments
- `isotherms`: `Vector` of `IsothermData` structs using the polymer in question with **only one kind of gas**. Temperature, pressure, density, and concentration must be provided.
- `bulk_parameters`: vector of the gas or vapor's characteristic parameters, in the format of `[p★_mpa, t★_k, ρ★_g_cm3]`.
- `polymer_parameters`: vector of the polymer's characteristic parameters, in the format of `[p★_mpa, t★_k, ρ★_g_cm3]`.
"""
function fit_kij(::NELF, isotherms::AbstractVector{<:IsothermData}, bulk_parameters::AbstractVector{<:Number}, polymer_parameters::AbstractVector{<:Number}; 
    interpolation_model=DualMode(), kij_fit_p_mpa=1e-4)
    
    interpolation_models = [fit_model(interpolation_model, isotherm) for isotherm in isotherms]
    bulk_model = SL(bulk_parameters...)
    densities = polymer_density.(isotherms) # get each isotherm's density in case the user accounted for polymers from different batches
    temperatures = temperature.(isotherms)
    polymer_phase_parameters = [[j for j in i] for i in zip(polymer_parameters, bulk_parameters)]  # hacky way to make it a vector. This is slow and dumb.
    function kij_error_function(kij)
        kij_mat = [0 kij[1]; kij[1] 0]
        ksw_vec = [0]
        polymer_model = SL(polymer_phase_parameters..., kij_mat)
        nelf_models = [NELFModel(bulk_model, polymer_model, densities[i], ksw_vec) for i in eachindex(isotherms)]
        pred_concs = [predict_concentration(nelf_models[i], temperatures[i], kij_fit_p_mpa, [1])[1] for i in eachindex(isotherms)]
        exp_concs = predict_concentration.(interpolation_models, kij_fit_p_mpa)
        return rss(exp_concs, pred_concs)
    end
    return Optim.minimizer(Optim.optimize(kij_error_function, [0.], BFGS()))
end


function fit_ksw(::NELF, isotherms::AbstractVector{<:IsothermData}, bulk_parameters::AbstractVector{<:Number}, polymer_parameters::AbstractVector{<:Number})
end

function fit_ksw(::NELF, isotherm::IsothermData, bulk_model, polymer_model)
    density = polymer_density(isotherm)
    temp = temperature(isotherm)
    pressures_mpa = partial_pressures(isotherm; component=1)
    concentrations_cc = concentration(isotherm; component=1)
    function error_function(ksw)
        nelf_model = NELFModel(bulk_model, polymer_model, density, ksw)
        pred_concs = [predict_concentration(nelf_model, temp, p, [1])[1] for p in pressures_mpa]
        return rss(concentrations_cc, pred_concs)
    end
    return Optim.minimizer(Optim.optimize(error_function, [0.], BFGS()))
end