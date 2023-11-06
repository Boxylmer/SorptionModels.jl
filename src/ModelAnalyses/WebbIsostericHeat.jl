
struct WebbIsostericHeatAnalysis <: AbstractIsostericHeatAnalysis
    isotherms
    sorption_models
    sampled_concentrations  # cc/cc
    temperatures # K
    inverse_temperature_vector  # ln(K)
    pressure_vectors # MPa
    ln_inv_vol_curves  # ln(1/(L/mol)) -> ln(mol/L)
    isosteric_heat_at_conc # J/mol
    isosteric_entropy_at_conc # j/(mol*K)
    pre_exponential_factors # (cc/cc) / MPa
end

"""
    WebbIsostericHeatAnalysis(isotherms::AbstractVector{<:IsothermData}, [eosmodel=missing]; 
        [model=DualMode()], 
        [num_points=25],
        [use_vant_hoff_constraints=false]
    ) 

Calculate the isosteric heat of sorption (``\\Delta{H}_{sorption}``) from a vector of isotherms as a function of concentration.
- `eosmodel`: A function that returns compressibility z when supplied a P (MPa) and T (K), i.e., z(P, T). When specified, will use this instead of assuming ideal behavior. 
- If the Dual Mode model is used, `use_vant_hoff_constraints` will constrain the dual mode fittings with respect to temperature. (see `VantHoffDualModeAnalysis`)
"""
function WebbIsostericHeatAnalysis(isotherms::AbstractVector{<:IsothermData}, eosmodel=missing; 
    model=DualMode(), num_points=25, 
    use_vant_hoff_constraints=false,
    )
    
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
    sampled_concentrations = LinRange(1e-3, max_concentration, num_points)

    # calculate the interpolated pressure curves and their logarithms (Base e)
    pressure_curves = [predict_pressure.(sorption_models, conc) for conc in sampled_concentrations]
    # first index -> pressure, second index -> temperature

    external_compressibilities = [
        [
            z(pressure_curves[p_index][t_index], temperatures[t_index]) 
            for t_index in eachindex(temperatures)
        ]
        for p_index in eachindex(pressure_curves)
    ]

    external_molar_volumes = [
        [
            external_compressibilities[p_index][t_index] * MembraneBase.R_MPA_L_K_MOL * temperatures[t_index] / pressure_curves[p_index][t_index] 
            for t_index in eachindex(temperatures)
        ]
        for p_index in eachindex(pressure_curves)
    ]

    for v in external_molar_volumes
        println(v)
    end
    # ln((cm3/mol / 1000cm3/L) / (L/mol)) = ln(L/mol / L/mol) -> no units

    # ln_inv_vol_curves = [(-MembraneBase.CC_PER_MOL_STP / 1000) ./ log.(mvs) for mvs in external_molar_volumes] 
    ln_inv_vol_curves = [-1 .* log.(mvs ./ (MembraneBase.CC_PER_MOL_STP / 1000)) for mvs in external_molar_volumes]

    # finally, actually get the slopes (linear fittings), and multiply by R to get the isosteric heats. 
    # If we're accounting for non-ideal terms, we multiply by Rz instead (where z is compressibility)
    isosteric_heat_at_conc = []
    isosteric_entropy_at_conc = []
    pre_exponential_factors = []


    for (i, ln_inv_vol_curve) in enumerate(ln_inv_vol_curves)
        slope, intercept =  fit_linear_data(inverse_temps, ln_inv_vol_curve)
        push!(isosteric_heat_at_conc, slope * MembraneBase.R_J_MOL_K)
        
        # Additionally, get the entropy of sorption (prototype, may be inaccurate)
        isosteric_entropy = MembraneBase.R_J_MOL_K * (log(sampled_concentrations[i]) - intercept)
        push!(isosteric_entropy_at_conc, isosteric_entropy)

        pre_exponential_factor = sampled_concentrations[i] * exp(-intercept)  # (cc/cc) / (MPa)
        push!(pre_exponential_factors, pre_exponential_factor)
    end

    return WebbIsostericHeatAnalysis(
        isotherms, sorption_models, sampled_concentrations,
        temperatures, inverse_temps, pressure_curves, ln_inv_vol_curves, 
        isosteric_heat_at_conc, isosteric_entropy_at_conc, pre_exponential_factors)
end