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

function bulk_phase_activity(model::NELFModel, temperature, pressure, bulk_phase_mole_fractions)
    μ = activity(
        model.bulk_model, 
        pressure, 
        temperature, 
        bulk_phase_mole_fractions)
    return μ
end

function predict_concentration(model::NELFModel, temperature, pressure, bulk_phase_mole_fractions; units=:cc)
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
- `isotherms`: `Vector` of `IsothermData` structs using the polymer in question. Temperature, pressure, and concentration must be provided.
- `bulk_phase_characteristic_params`: Vector of pure characteristic parameter vectors following the same order as the isotherms. 
  - E.g., for Sanchez Lacombe and two input isotherms, `bulk_phase_characteristic_params = [[p★_1, t★_1, ρ★_1, mw_1], [p★_2, t★_2, ρ★_2, mw_2]]`
"""
fit_model(::NELF, ::SL, isotherms::AbstractVector{<:IsothermData}, bulk_phase_characteristic_params; 
    polymer_molecular_weight=100000)
    
    # this function uses SL, which needs 4 params per component, one of which is already specified (MW)
    
    dualmode_models = fit_model.(DualMode(), isotherms)
    
    bulk_phase_models = [SL(params...) for params in bulk_phase_characteristic_params]
    default_ksw_vec = [0]
    densities = polymer_density.(isotherms) # get each isotherm's density in case the user accounted for polymers from different batches
    infinite_dilution_pressure = eps() # ???
    
    function error_function(char_param_vec)
        polymer_phase_models = [SL(zip(char_param_vec, params)...) for params in bulk_phase_characteristic_params]
        nelf_models = [NELFModel(bulk_phase_models[i], polymer_phase_models[i], densities[i], default_ksw_vec)]
        predict_concentration(nelf_models, )
    end
    # work in progress

end