struct NELF end


"""
    Requires that in the kijmatrix, the polymer is the first index. The remaining indexes must match `penetrants`.
"""

struct NELFModel{BMT, POLYMT, PDT} <: SorptionModel
    bulk_model::BMT                 # EOS
    polymer_model::POLYMT           # EOS
    polymer_dry_density::PDT        # number
    # ksw_values::KSWT                # vector of values
end

function predict_concentration(model::NELFModel, temperature::Number, pressure::Number, bulk_phase_mole_fractions; ksw=nothing, units=:cc)
    minimum_val = 100 * eps()
    error_target = _make_nelf_model_mass_fraction_target(model, temperature, pressure, bulk_phase_mole_fractions; ksw, minimum_val)

    penetrant_mass_fraction_initial_guesses = ones(length(bulk_phase_mole_fractions)) * eps()
    lower = zeros(length(penetrant_mass_fraction_initial_guesses))
    upper = ones(length(penetrant_mass_fraction_initial_guesses)) .- eps()
    res = Optim.optimize(
        error_target, 
        lower, upper,
        penetrant_mass_fraction_initial_guesses, 
        # Fminbox(LBFGS()),
        Fminbox(NelderMead()),
        # Newton(),
        # SAMIN(),
        Optim.Options(
        allow_f_increases = false,
        # x_tol = 1e-9,
        # g_tol = 1e-7
    ); autodiff=:forward)

    penetrant_mass_fractions = Optim.minimizer(res)
    polymer_phase_mass_fractions = vcat(1 - sum(penetrant_mass_fractions), penetrant_mass_fractions)

    if units==:frac
        return polymer_phase_mass_fractions[2:end]
    elseif units==:g
        concs_g_g = polymer_phase_mass_fractions_to_gpen_per_gpol(polymer_phase_mass_fractions)
        return concs_g_g
    elseif units==:cc
        concs_cc_cc = polymer_phase_mass_fractions_to_ccpen_per_ccpol(polymer_phase_mass_fractions, model.polymer_dry_density, molecular_weight(model.bulk_model))
        return concs_cc_cc
    end
end
predict_concentration(model::NELFModel, temperature::Number, pressure::Number; kwargs...) = predict_concentration(model, temperature, pressure, [1]; kwargs...)

function _make_nelf_model_mass_fraction_target(model::NELFModel, temperature::Number, pressure::Number, bulk_phase_mole_fractions; ksw=nothing, minimum_val=100*eps())
    if isnothing(ksw)
        ksw = zeros(length(bulk_phase_mole_fractions))
    end

    pressure = pressure < minimum_val ? minimum_val : pressure    
    bulk_phase_mole_fractions = ((i) -> (i < minimum_val ? minimum_val : i)).(bulk_phase_mole_fractions)  # Courtesy of Clementine (Julia Discord)

    target_activities = chemical_potential(model.bulk_model, pressure, temperature, bulk_phase_mole_fractions)
    normalizer = (rss(zeros(length(bulk_phase_mole_fractions)), target_activities))

    function error_target(penetrant_mass_fractions)
        polymer_mass_fraction = 1 - sum(penetrant_mass_fractions)
        if polymer_mass_fraction <= 0 || polymer_mass_fraction > 1
            @warn "Polymer mass fraction was not valid, returning very large error value: " * string(1 - sum(penetrant_mass_fractions))
            return log1p(normalizer)
            return 1e100
        end
        polymer_phase_mass_fractions = vcat(polymer_mass_fraction, penetrant_mass_fractions)
        polymer_phase_density_after_swelling = calculate_polymer_phase_density(model, pressure, bulk_phase_mole_fractions, polymer_phase_mass_fractions, ksw)
        polymer_phase_density_upper_bound = density_upper_bound(model.polymer_model, polymer_phase_mass_fractions)
        if polymer_phase_density_after_swelling > polymer_phase_density_upper_bound 
            return log1p(normalizer)
            return 1e100
        end

        polymer_phase_activities = ρTω_chemical_potential(
            model.polymer_model, 
            polymer_phase_density_after_swelling, 
            temperature, 
            polymer_phase_mass_fractions)
        residual_squared = log1p((rss(target_activities, polymer_phase_activities[2:end])))
        return residual_squared
    end
    return error_target
end


function calculate_swelled_polymer_density(model::NELFModel, penetrant_partial_pressures::AbstractVector{<:Number}, ksw_values::AbstractVector{<:Number})
    return model.polymer_dry_density * (1 - sum(ksw_values .* penetrant_partial_pressures))
end

function calculate_swelled_polymer_density(model::NELFModel, pressure::Number, penetrant_mole_fractions::AbstractVector{<:Number}, ksw_values::AbstractVector{<:Number})
    partial_pressures = penetrant_mole_fractions .* pressure
    return calculate_swelled_polymer_density(model, partial_pressures, ksw_values)
end

function calculate_polymer_phase_density(swelled_density::Number, polymer_phase_mass_fractions::AbstractVector{<:Number})
    return swelled_density / polymer_phase_mass_fractions[1]
end

function calculate_polymer_phase_density(model::NELFModel, pressure::Number, bulk_phase_mole_fractions::AbstractVector{<:Number}, polymer_phase_mass_fractions::AbstractVector{<:Number}, ksw_values::AbstractVector{<:Number})
    polymer_density_after_swelling = calculate_swelled_polymer_density(model, pressure, bulk_phase_mole_fractions, ksw_values)
    return calculate_polymer_phase_density(polymer_density_after_swelling, polymer_phase_mass_fractions)
end

"""
    infinite_dilution_solubility(model::NELFModel, temperature::Number)
Currenlty only supported for Sanchez Lacombe based models, get infinite dilution solubility in **((CC/CC) / MPa)**
"""
function infinite_dilution_solubility(model::NELFModel, temperature::Number)
    if typeof(model.polymer_model) <: MembraneEOS.SanchezLacombeModel
        t_st = 273.15 # K
        p_st = 0.1  # MPa
        comps = model.polymer_model.components
        polymer = comps[1]
        penetrant = comps[2]
        p★_pol, t★_pol, ρ★_pol = MembraneEOS.characteristic_pressure(polymer), MembraneEOS.characteristic_temperature(polymer), MembraneEOS.characteristic_density(polymer)
        p★_pen, t★_pen, ρ★_pen, mw_pen = MembraneEOS.characteristic_pressure(penetrant), MembraneEOS.characteristic_temperature(penetrant), MembraneEOS.characteristic_density(penetrant), MembraneEOS.molecular_weight(penetrant)
        p★12 = (sqrt(p★_pol) - sqrt(p★_pen))^2
        ρ_pol = model.polymer_dry_density
        R = MembraneBase.R_MPA_L_K_MOL * 1000  # --> MPa * cm3 / molK
        
        coeff = t_st / (temperature * p_st)
        exp_term_1 = mw_pen * p★_pen / (ρ★_pen * R * t★_pen) * (1 + ((t★_pen * p★_pol)/(t★_pol * p★_pen) - 1) * (ρ★_pol / ρ_pol)) * log1p(-ρ_pol/ρ★_pol)
        exp_term_2 = (t★_pen * p★_pol)/(t★_pol * p★_pen) - 1
        exp_term_3 = (ρ_pol * t★_pen)/(ρ★_pol * p★_pen * temperature) * (p★_pen + p★_pol - p★12)

        s_inf = coeff * exp(exp_term_1 + exp_term_2 + exp_term_3)
        return s_inf
        # if length(model.polymer_model.components) > 2
        #     throw(ErrorException("NELF model must be a pure component to "))
        # end
    else
        throw(ErrorException("Only NELF models using Sanchez Lacombe support infinite dilution at this time"))
    end


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
    polymer_molecular_weight=100000, eos=MembraneEOS.SanchezLacombe())
    
    # this function uses SL, which needs 4 params per component, one of which is already specified (MW)

    infinite_dilution_pressure = 1e-5 # ???

    error_function = _make_nelf_model_parameter_target(isotherms, bulk_phase_characteristic_params, infinite_dilution_pressure, polymer_molecular_weight)
    
    densities = polymer_density.(isotherms)
    density_lower_bound = maximum(densities)
    
    lower = [0., 0., density_lower_bound]
    upper = [3000, 3000, 3.]
    res = Optim.optimize(
        error_function, lower, upper, 
        [500, 500, upper[3] - 100 * eps()], 
        # Fminbox(LBFGS(; m=60, linesearch = Optim.LineSearches.BackTracking())), 
        Fminbox(LBFGS()), 
        # Fminbox(NelderMead()),
        Optim.Options(; allow_f_increases = true))
    # res = Optim.optimize(
    #     error_function, lower, upper, 
    #     [100, 100, density_lower_bound * 1.2], 
    #     SAMIN(; rt = 0.04), 
    #     Optim.Options(iterations=10^6))

    return [Optim.minimizer(res)..., polymer_molecular_weight]
    # work in progress
end
function _make_nelf_model_parameter_target(isotherms, bulk_phase_characteristic_params, infinite_dilution_pressure, polymer_molecular_weight=100000)
    bulk_phase_models = [SL(params...) for params in bulk_phase_characteristic_params]
    dualmode_models = [fit_model(DualMode(), isotherm) for isotherm in isotherms]
    densities = polymer_density.(isotherms) # get each isotherm's density in case the user accounted for polymers from different batches
    temperatures = temperature.(isotherms)

    function error_function(char_param_vec)
        given_sol = zeros(length(isotherms))
        pred_sol = zeros(length(isotherms)) 
        for i in eachindex(isotherms)
            char_pressures = [char_param_vec[1], bulk_phase_characteristic_params[i][1]]
            char_temperatures = [char_param_vec[2], bulk_phase_characteristic_params[i][2]]
            char_densities = [char_param_vec[3], bulk_phase_characteristic_params[i][3]]
            molecular_weights = [polymer_molecular_weight, bulk_phase_characteristic_params[i][4]]

            polymer_phase_model = SL(char_pressures, char_temperatures, char_densities, molecular_weights)
            nelf_model = NELFModel(bulk_phase_models[i], polymer_phase_model, densities[i])
            pred_sol[i] = predict_concentration(nelf_model, temperatures[i], infinite_dilution_pressure, [1]; ksw=[0])[1] / infinite_dilution_pressure
            given_sol[i] = infinite_dilution_solubility(dualmode_models[i]::DualModeModel)
        end
        err = log1p(rss(given_sol, pred_sol))
        # @show char_param_vec, err
        return err
    end
    return error_function
end

function _make_nelf_model_parameter_target_2(isotherms, bulk_phase_characteristic_params, infinite_dilution_pressure, polymer_molecular_weight=100000)
    bulk_phase_models = [SL(params...) for params in bulk_phase_characteristic_params]
    dualmode_models = [fit_model(DualMode(), isotherm) for isotherm in isotherms]
    densities = polymer_density.(isotherms) # get each isotherm's density in case the user accounted for polymers from different batches
    temperatures = temperature.(isotherms)

    function error_function(char_param_vec)
        given_sol_inf = zeros(length(isotherms))  # reused
        pred_sol_inf = zeros(length(isotherms))   # reused

        pred_errs = 0

        for i in eachindex(isotherms)
            char_pressures = [char_param_vec[1], bulk_phase_characteristic_params[i][1]]
            char_temperatures = [char_param_vec[2], bulk_phase_characteristic_params[i][2]]
            char_densities = [char_param_vec[3], bulk_phase_characteristic_params[i][3]]
            molecular_weights = [polymer_molecular_weight, bulk_phase_characteristic_params[i][4]]

            polymer_phase_model = SL(char_pressures, char_temperatures, char_densities, molecular_weights)
            nelf_model = NELFModel(bulk_phase_models[i], polymer_phase_model, densities[i])

            pressures = partial_pressures(isotherms[i]; component=1)
            if pressures[1] == 0; pressures = partial_pressures(isotherms[i]; component=1)[2:end]; end

            pred_sols = [predict_concentration(nelf_model, temperatures[i], p, [1]; ksw=[0])[1] / p for p in pressures]
            given_sols = [predict_concentration(dualmode_models[i]::DualModeModel, p) / p for p in pressures]
            pred_errs += log(rss(given_sols, pred_sols))

            pred_sol_inf[i] = predict_concentration(nelf_model, temperatures[i], infinite_dilution_pressure, [1]; ksw=[0])[1] / infinite_dilution_pressure
            given_sol_inf[i] = infinite_dilution_solubility(dualmode_models[i]::DualModeModel)
        end
        err_inf = log(rss(given_sol_inf, pred_sol_inf))
        
        # @show char_param_vec, err
        return pred_errs
        return err_inf + pred_errs
    end
    return error_function
end

function _make_nelf_model_parameter_target_3(isotherms, bulk_phase_characteristic_params, infinite_dilution_pressure, polymer_molecular_weight=100000)
    bulk_phase_models = [SL(params...) for params in bulk_phase_characteristic_params]
    dualmode_models = [fit_model(DualMode(), isotherm) for isotherm in isotherms]
    densities = polymer_density.(isotherms) # get each isotherm's density in case the user accounted for polymers from different batches
    temperatures = temperature.(isotherms)

    function error_function(char_param_vec)
        given_sol = zeros(length(isotherms))  # reused
        pred_sol = zeros(length(isotherms))   # reused
        for i in eachindex(isotherms)
            char_pressures = [char_param_vec[1], bulk_phase_characteristic_params[i][1]]
            char_temperatures = [char_param_vec[2], bulk_phase_characteristic_params[i][2]]
            char_densities = [char_param_vec[3], bulk_phase_characteristic_params[i][3]]
            molecular_weights = [NaN, bulk_phase_characteristic_params[i][4]]

            polymer_phase_model = SL(char_pressures, char_temperatures, char_densities, molecular_weights)
            nelf_model = NELFModel(bulk_phase_models[i], polymer_phase_model, densities[i])
            pred_sol[i] = (infinite_dilution_solubility(nelf_model, temperatures[i]))
            given_sol[i] = (infinite_dilution_solubility(dualmode_models[i]::DualModeModel))
        end
        err = log1p(rss(given_sol, pred_sol))
        # @show char_param_vec, err
        return err
    end
    return error_function
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
        nelf_models = [NELFModel(bulk_model, polymer_model, densities[i]) for i in eachindex(isotherms)]
        pred_concs = [predict_concentration(nelf_models[i], temperatures[i], kij_fit_p_mpa, [1]; ksw=ksw_vec)[1] for i in eachindex(isotherms)]
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
    nelf_model = NELFModel(bulk_model, polymer_model, density)
    function error_function(ksw)
        pred_concs = [predict_concentration(nelf_model, temp, p, [1]; ksw=ksw)[1] for p in pressures_mpa]
        @show ksw
        return rss(concentrations_cc, pred_concs)
    end
    # lower = [0.]
    # upper = [Inf]
    init = [0.0]
    @show res = Optim.optimize(error_function, init, BFGS())
    @show minimizer = Optim.minimizer(res)
    return min(0.0, minimizer)
end