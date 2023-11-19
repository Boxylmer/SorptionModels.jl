abstract type AbstractIsostericHeatAnalysis end

struct IsostericHeatAnalysis <: AbstractIsostericHeatAnalysis
    isotherms
    sorption_models
    sampled_concentrations  # cc/cc
    temperatures # K
    inverse_temperature_vector  # ln(K)
    pressure_vectors # MPa
    ln_pressure_vectors  # ln(mpa)
    isosteric_heat_at_conc # J/mol
    isosteric_entropy_at_conc # j/(mol*K)
    pre_exponential_factors # (cc/cc) / MPa
    z_values
end

"""
    IsostericHeatAnalysis(isotherms::AbstractVector{<:IsothermData}, [eosmodel=missing]; 
        [model=DualMode()], 
        [num_points=25],
        [use_vant_hoff_constraints=false],
        [gab_pressure_conversion_funcs=missing], 
        [gab_activity_conversion_funcs=missing],
        [linear_regression_error=true]) 

Calculate the isosteric heat of sorption (``\\Delta{H}_{sorption}``) from a vector of isotherms as a function of concentration.
- `eosmodel`: A function that returns compressibility z when supplied a P (MPa) and T (K), i.e., z(P, T). When specified, will use this instead of assuming ideal behavior. 
- If the Dual Mode model is used, `use_vant_hoff_constraints` will constrain the dual mode fittings with respect to temperature. (see `VantHoffDualModeAnalysis`)
- If the GAB model is used, `gab_pressure_conversion_funcs` (converting pressure to activity) and `gab_activity_conversion_funcs` (convert activity to pressure) will need to be defined in the same order as the isotherms.
- if `linear_regression_error` is true, each concentration will have uncertainties (Measurements.jl) calculated using the standard error of linear regression. If any inputs have uncertainty, downstream calculations' uncertainties will be used to weight the regression.  
"""

function IsostericHeatAnalysis(isotherms::AbstractVector{<:IsothermData}, eosmodel=missing; 
    model=DualMode(), num_points=25, 
    use_vant_hoff_constraints=false,
    gab_pressure_conversion_funcs=missing, 
    gab_activity_conversion_funcs=missing,
    linear_regression_error=true)
    
    if ismissing(eosmodel)
        z(p, t) = 1
    else
        z = eosmodel
    end

    # first fit some models to each isotherm so that we can accurately interpolate
    if use_vant_hoff_constraints  && typeof(model) <: DualMode
        sorption_models = Vector(
            strip_measurement_to_value(
                VantHoffDualModeAnalysis(isotherms; use_fugacity=false).final_models
            )
        )
    elseif typeof(model) <: GAB
        if ismissing(gab_pressure_conversion_funcs)
            throw(ErrorException("The GAB model supplied has no way of converting pressures back to activities."))
        end
        if ismissing(gab_activity_conversion_funcs)
            throw(ErrorException("The GAB model supplied has no way of converting activities to pressures."))
        end
        sorption_models = [fit_model(
            model, isotherms[idx]; 
            uncertainty_method=nothing, 
            pressure_conversion_function=gab_pressure_conversion_funcs[idx],
            activity_conversion_function=gab_activity_conversion_funcs[idx]) 
            for idx in 1:length(isotherms)]

    else
        sorption_models = [fit_model(model, isotherm; uncertainty_method=nothing) for isotherm in isotherms]
    end
    # calculate the x axis of the isosteric heat slope (inverse temperature)
    temperatures = temperature.(isotherms)
    inverse_temps = 1 ./ temperatures
    
    # Get the concentrations in each isotherm to find which one has the lowest maximum concentration, so that we know how far we can interpolate
    isotherm_concentration_vectors = [concentration(isotherm; component=1) for isotherm in isotherms]
    max_concentration = minimum([isotherm_concentration_vector[end] for isotherm_concentration_vector in isotherm_concentration_vectors])
    # and pick a set of concentrations to sample at below that maximum
    sampled_concentrations = LinRange(1e-1, max_concentration, num_points)
    # calculate the interpolated pressure curves and their logarithms (Base e)
    pressure_curves = [predict_pressure.(sorption_models, conc) for conc in sampled_concentrations]
    ln_pressure_curves = [log.(pressure_curve) for pressure_curve in pressure_curves]
    
    
    
    
    all_z_values = [[z(pressure_curves[p_idx][t_idx], temperatures[t_idx]) for t_idx in eachindex(temperatures)] for p_idx in eachindex(pressure_curves)]
    avg_z_values = sum.(all_z_values) ./ length.(all_z_values)
    
    # finally, actually get the slopes (linear fittings), and multiply by R to get the isosteric heats. 
    # If we're accounting for non-ideal terms, we multiply by Rz instead (where z is compressibility)
    if linear_regression_error
        VALTYPE = Measurement
    else
        VALTYPE = eltype(promote(all_z_values[1][1], sampled_concentrations[1]))
    end

    isosteric_heat_at_conc = VALTYPE[]
    isosteric_entropy_at_conc = VALTYPE[]
    pre_exponential_factors = VALTYPE[]
    
    for (i, ln_pressure_curve) in enumerate(ln_pressure_curves)
        slope, intercept = fit_linear_data(inverse_temps, ln_pressure_curve)
        if !linear_regression_error
            slope = strip_measurement_to_value(slope)
            intercept = strip_measurement_to_value(intercept)
        end
        push!(isosteric_heat_at_conc, slope * avg_z_values[i] * MembraneBase.R_J_MOL_K)
        
        # Additionally, get the entropy of sorption (prototype, may be inaccurate)
        # TODO does this need z as well?
        isosteric_entropy = MembraneBase.R_J_MOL_K * (log(sampled_concentrations[i]) - intercept)
        push!(isosteric_entropy_at_conc, isosteric_entropy)

        pre_exponential_factor = sampled_concentrations[i] * exp(-intercept)  # (cc/cc) / (MPa)
        push!(pre_exponential_factors, pre_exponential_factor)
    end

    return IsostericHeatAnalysis(
        isotherms, sorption_models, sampled_concentrations,
        temperatures, inverse_temps, pressure_curves, ln_pressure_curves, 
        isosteric_heat_at_conc, isosteric_entropy_at_conc, pre_exponential_factors,
        all_z_values)
end