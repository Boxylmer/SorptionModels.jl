
struct IsostericHeatAnalysis
    isotherms
    sorption_models
    sampled_concentrations  # cc/cc
    temperatures # K
    inverse_temperature_vector  # ln(K)
    pressure_vectors # MPa
    ln_pressure_vectors  # ln(mpa)
    isosteric_heat_at_conc # J/mol
end

"""
    (isotherms::AbstractVector{<:IsothermData}; model=DualMode(), num_points=25) 

Calculate the isosteric heat of sorption (``\\Delta{H}_{sorption}``) from a vector of isotherms as a function of concentration.
"""
function IsostericHeatAnalysis(isotherms::AbstractVector{<:IsothermData}; model=DualMode(), num_points=25, use_vant_hoff_constraints=false)
    # first fit some models to each isotherm so that we can accurately interpolate
    if use_vant_hoff_constraints
        sorption_models = Vector(
            strip_measurement_to_value(
                VantHoffDualModeAnalysis(isotherms; use_fugacity=false).final_models
            )
        )

    else
        sorption_models = [fit_model(model, isotherm; uncertainty_method=:JackKnife) for isotherm in isotherms]
    end
    # calculate the x axis of the isosteric heat slope (inverse temperature)
    temperatures = temperature.(isotherms)
    inverse_temps = 1 ./ temperatures
    
    # Get the concentrations in each isotherm to find which one has the lowest maximum concentration, so that we know how far we can interpolate
    isotherm_concentration_vectors = [concentration(isotherm; component=1) for isotherm in isotherms]
    max_concentration = minimum([isotherm_concentration_vector[end] for isotherm_concentration_vector in isotherm_concentration_vectors])
    # and pick a set of concentrations to sample at below that maximum
    sampled_concentrations = LinRange(1e-3, max_concentration, num_points)
    # calculate the interpolated pressure curves and their logarithms (Base e)
    pressure_curves = [predict_pressure.(sorption_models, conc) for conc in sampled_concentrations]
    ln_pressure_curves = [log.(pressure_curve) for pressure_curve in pressure_curves]

    # finally, actually get the slopes (linear fittings), and multiply by R to get the isosteric heats. 
    # If we're accounting for non-ideal terms, we multiply by Rz instead (where z is compressibility)
    isosteric_heat_at_conc = []
    for ln_pressure_curve in ln_pressure_curves
        slope, _ =  fit_linear_data(inverse_temps, ln_pressure_curve)
        push!(isosteric_heat_at_conc, slope * MembraneBase.R_J_MOL_K)
    end

    return IsostericHeatAnalysis(
        isotherms, sorption_models, sampled_concentrations,
        temperatures, inverse_temps, pressure_curves, ln_pressure_curves, 
        isosteric_heat_at_conc)
end