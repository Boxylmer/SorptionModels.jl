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
    MolarVolumeAnalysis(model::SorptionModel, pressures_mpa, frac_dilations, [uncertainty_method=:Hessian], [n_interp=30], [n_params=3])

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
function MolarVolumeAnalysis(model::SorptionModel, pressures_mpa::AbstractVector{<:Number}, frac_dilations::AbstractVector{<:Number}, 
    isothermal_compressability=0; uncertainty_method=:Hessian, n_interp=30, n_params = 3)
    
    concentrations = predict_concentration(model, pressures_mpa)

    continuous_pressure_curve(c_ccpercc) = predict_pressure(model, c_ccpercc)

    dp_dc = ForwardDiff.derivative.(continuous_pressure_curve, concentrations) # mpa / cc/cc

    dilation_function_params = find_dilation_function_params(pressures_mpa, frac_dilations, uncertainty_method; start = ones(n_params))

    continuous_dilations = dilation_empirical_function.(pressures_mpa, dilation_function_params...)
    continuous_dilation_derivatives = ForwardDiff.derivative.(x -> dilation_empirical_function(x, dilation_function_params...), pressures_mpa)

    diagnostic_pressures_mpa = collect(range(minimum(pressures_mpa), maximum(pressures_mpa), n_interp))
    diagnostic_frac_dilations = dilation_empirical_function.(diagnostic_pressures_mpa, dilation_function_params...)
    diagnostic_concentrations = predict_concentration(model, diagnostic_pressures_mpa)
    diagnostic_dp_dc = ForwardDiff.derivative.(continuous_pressure_curve, diagnostic_concentrations)
    diagnostic_dilation_derivatives = ForwardDiff.derivative.(x -> dilation_empirical_function(x, dilation_function_params...), diagnostic_pressures_mpa)
    diagnostic_volumes = (diagnostic_dilation_derivatives .+ isothermal_compressability) .* diagnostic_dp_dc .* MembraneBase.CC_PER_MOL_STP # (1/MPa) * (Mpa / (cc/cc)) * (cc/mol) = cm3/mol

    volumes = (continuous_dilation_derivatives .+ isothermal_compressability) .* dp_dc .* MembraneBase.CC_PER_MOL_STP # (1/MPa) * (Mpa / (cc/cc)) * (cc/mol) = cm3/mol
    return MolarVolumeAnalysis(pressures_mpa, concentrations, dp_dc, frac_dilations, continuous_dilations, continuous_dilation_derivatives, volumes,
        diagnostic_pressures_mpa, diagnostic_frac_dilations, diagnostic_concentrations, diagnostic_dp_dc, diagnostic_dilation_derivatives, diagnostic_volumes)
end

function dilation_empirical_function(p_mpa, a, b, c, d = 0)
    res = a * p_mpa /(1 + b*p_mpa) + c * p_mpa + d * p_mpa^2
    return res
end


function dilation_empirical_function_sqerr(X, Y, params)
    sum((dilation_empirical_function.(X, params...) .- Y).^2)
end

function find_dilation_function_params(pressures_mpa, frac_dilations, uncertainty_method=nothing; start = [1., 1, 1])
    obj = x -> dilation_empirical_function_sqerr(
        strip_measurement_to_value(pressures_mpa), 
        strip_measurement_to_value(frac_dilations), 
        x
    )
    n = length(pressures_mpa)
    res = Optim.optimize(obj, start, BFGS()).minimizer

    if uncertainty_method == :Hessian
        model_uncertainty = rss_minimizer_standard_errors(obj, res, n)
        return res .± model_uncertainty

    elseif uncertainty_method == :JackKnife
        data = collect(zip(pressures_mpa, frac_dilations))
        fit_function(data) = find_dilation_function_params(collect.(collect(zip(data...)))...; start = res)
        model_uncertainty = jackknife_uncertainty(fit_function, data)
        return res .± model_uncertainty

    else
        return res
    end
end


