struct MolarVolumeAnalysis
    pressures_mpa
    concentrations_cc_cc
    dp_dc
    frac_dilations
    continuous_dilations
    dfracional_dilation_dp
    partial_molar_volumes_cm3_mol

    diagnostic_pressures_mpa
    diagnostic_frac_dilations
    diagnostic_concentrations
    diagnostic_dp_dc
    diagnostic_dilation_derivatives
    diagnostic_volumes
end

"""
    MolarVolumeAnalysis(model::SorptionModel, pressures_mpa, frac_dilations, [uncertainty_method=nothing], [n_interp=30], [n_params=3])

Apply the calculation discussed in

R. Raharjo, B. Freeman, E. Sanders, Pure and mixed gas CH4 and n-C4H10 sorption and dilation in poly(dimethylsiloxane), Journal of Membrane Science. 292 (2007) 45–61. https://doi.org/10.1016/j.memsci.2007.01.012.

to calculate the partial molar volume of a component in a polymer phase.

- The model in question should take true pressures and not fugacities. 
- The isothermal compressibility factor (units of MPa^-1) is neglected by default. It is used for calculating the change in volume due to external pressure and can generally be neglected for condensible gasses, low pressure liquids, and vapors. For permanent gasses and high pressure liquids, ensure this can be neglected or specify it's value.
- A modified dual mode model is used to approximate the derivative of the dilation data. If you're uncertainty is very high and you're using :Hessian uncertainties, it's likely you're overfitting and should use :JackKnife instead, as it will more adequately represent the uncertainty in your fit.  

- uncertainty_method
- n_interp will determine the number of values to interpolate to
- n_params will specify the number of fitted parameters to use when fitting a continuous dilation function. Currently only 3 and 4 are supported. Start with 3, and if the dilation data is non-monotonic or you're dealing with some difficult data, move to 4. Be aware this will increase the degrees of freedom in the fitting and thus the uncertainty associated with it (see `uncertainty_method`).
"""
function MolarVolumeAnalysis(sorptionmodel::SorptionModel, pressures_mpa::AbstractVector{<:Number}, frac_dilations::AbstractVector{<:Number}, 
    isothermal_compressability=0; uncertainty_method=nothing, n_interp=30, n_params = 3, modeltype=EmpiricalDilation())
    
    concentrations = predict_concentration(sorptionmodel, pressures_mpa)

    continuous_pressure_curve(c_ccpercc) = predict_pressure(sorptionmodel, c_ccpercc)

    dp_dc = ForwardDiff.derivative.(continuous_pressure_curve, concentrations) # mpa / cc/cc

    # prev 
    # dilation_function_params = find_dilation_function_params(pressures_mpa, frac_dilations, uncertainty_method; start = ones(n_params))
    # continuous_dilations = dilation_empirical_function.(pressures_mpa, dilation_function_params...)
    # continuous_dilation_derivatives = ForwardDiff.derivative.(x -> dilation_empirical_function(x, dilation_function_params...), pressures_mpa)
    # diagnostic_pressures_mpa = collect(range(minimum(pressures_mpa), maximum(pressures_mpa), n_interp))
    # diagnostic_frac_dilations = dilation_empirical_function.(diagnostic_pressures_mpa, dilation_function_params...)
    
    # new 
    if modeltype == EmpiricalDilation()
        dilation_model = fit_model(modeltype, pressures_mpa, frac_dilations, uncertainty_method; n_params)
    elseif modeltype == DualModeDilation()
        dilation_model = fit_model(modeltype, pressures_mpa, frac_dilations, sorptionmodel, uncertainty_method)
    else
        throw(ArgumentError("Dilation model type not supported."))
    end
    
    continuous_dilations = predict_dilation(dilation_model, pressures_mpa)
    continuous_dilation_derivatives = predict_dilation_derivative(dilation_model, continuous_dilations)
    # ForwardDiff.derivative.(x -> predict_dilation(dilation_model, x), pressures_mpa)
    diagnostic_pressures_mpa = collect(range(minimum(pressures_mpa), maximum(pressures_mpa), n_interp))
    diagnostic_frac_dilations = predict_dilation(dilation_model, diagnostic_pressures_mpa)


    diagnostic_concentrations = predict_concentration(sorptionmodel, diagnostic_pressures_mpa)
    diagnostic_dp_dc = ForwardDiff.derivative.(continuous_pressure_curve, diagnostic_concentrations)


    # prev
    # diagnostic_dilation_derivatives = ForwardDiff.derivative.(x -> dilation_empirical_function(x, dilation_function_params...), diagnostic_pressures_mpa)
    # diagnostic_volumes = (diagnostic_dilation_derivatives .+ isothermal_compressability) .* diagnostic_dp_dc .* MembraneBase.CC_PER_MOL_STP # (1/MPa) * (Mpa / (cc/cc)) * (cc/mol) = cm3/mol

    # new
    diagnostic_dilation_derivatives = predict_dilation_derivative(dilation_model, diagnostic_pressures_mpa)
    diagnostic_volumes = (diagnostic_dilation_derivatives .+ isothermal_compressability) .* diagnostic_dp_dc .* MembraneBase.CC_PER_MOL_STP # (1/MPa) * (Mpa / (cc/cc)) * (cc/mol) = cm3/mol


    volumes = (continuous_dilation_derivatives .+ isothermal_compressability) .* dp_dc .* MembraneBase.CC_PER_MOL_STP # (1/MPa) * (Mpa / (cc/cc)) * (cc/mol) = cm3/mol
    return MolarVolumeAnalysis(pressures_mpa, concentrations, dp_dc, frac_dilations, continuous_dilations, continuous_dilation_derivatives, volumes,
        diagnostic_pressures_mpa, diagnostic_frac_dilations, diagnostic_concentrations, diagnostic_dp_dc, diagnostic_dilation_derivatives, diagnostic_volumes)
end
