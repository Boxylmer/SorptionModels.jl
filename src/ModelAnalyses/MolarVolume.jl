struct MolarVolumeAnalysis{PT, CT, DPDCT, FDT, CDT, DFDPT, PMVT}
    pressures_mpa::PT
    concentrations_cc_cc::CT
    dp_dc::DPDCT
    frac_dilations::FDT
    continuous_dilations::CDT
    dfracional_dilation_dp::DFDPT
    partial_molar_volumes_cm3_mol::PMVT

    diagnostic_pressures_mpa::PT
    diagnostic_frac_dilations::CDT
end

"""
    MolarVolumeAnalysis(model::SorptionModel, pressures_mpa, frac_dilations)

Apply the calculation discussed in

R. Raharjo, B. Freeman, E. Sanders, Pure and mixed gas CH4 and n-C4H10 sorption and dilation in poly(dimethylsiloxane), Journal of Membrane Science. 292 (2007) 45–61. https://doi.org/10.1016/j.memsci.2007.01.012.

to calculate the partial molar volume of a component in a polymer phase.

- The model in question should take true pressures and not fugacities. 
- The isothermal compressibility factor (units of MPa^-1) is neglected by default. It is used for calculating the change in volume due to external pressure and can generally be neglected for condensible gasses, low pressure liquids, and vapors. For permanent gasses and high pressure liquids, ensure this can be neglected or specify it's value.
- Polynomials are currently used to approximate the derivative of the dilation data, `poly_fit` specifies the degree of the polynomial to be used. 
"""
function MolarVolumeAnalysis(model::SorptionModel, pressures_mpa::AbstractVector{<:Number}, frac_dilations::AbstractVector{<:Number}, 
    isothermal_compressability=0, poly_fit=4, n_interp=30)
    
    concentrations = predict_concentration(model, pressures_mpa)

    continuous_pressure_curve(c_ccpercc) = predict_pressure(model, c_ccpercc)

    dp_dc = ForwardDiff.derivative.(continuous_pressure_curve, concentrations) # mpa / cc/cc
    
    # continuous_dilation_curve = fit(pressures_mpa, frac_dilations, poly_fit)
    # continuous_dilation_curve_derivative = derivative(continuous_dilation_curve, 1)
    # continuous_dilations = continuous_dilation_curve.(pressures_mpa)
    # continuous_dilation_derivatives = continuous_dilation_curve_derivative.(pressures_mpa)

    dilation_function_params = find_dilation_function_params(pressures_mpa, frac_dilations, :Hessian)
    continuous_dilations = dilation_empirical_function.(pressures_mpa, dilation_function_params...)
    continuous_dilation_derivatives = ForwardDiff.derivative.(x -> dilation_empirical_function(x, dilation_function_params...), pressures_mpa)

    diagnostic_pressures_mpa = collect(range(minimum(pressures_mpa), maximum(pressures_mpa), n_interp))
    diagnostic_frac_dilations = dilation_empirical_function.(diagnostic_pressures_mpa, dilation_function_params...)

    # dfracdil_dp = estimate_slope_by_adjacent_points(pressures_mpa, frac_dilations) # 1 / mpa

    volumes = (continuous_dilation_derivatives .+ isothermal_compressability) .* dp_dc .* MembraneBase.CC_PER_MOL_STP # (1/MPa) * (Mpa / (cc/cc)) * (cc/mol) = cm3/mol
    return MolarVolumeAnalysis(pressures_mpa, concentrations, dp_dc, frac_dilations, continuous_dilations, continuous_dilation_derivatives, volumes,
        diagnostic_pressures_mpa, diagnostic_frac_dilations)
end

# function dilation_empirical_function(p_mpa, a, b, c)
#     res = a * b * c * p_mpa / ((1 - b * p_mpa) * (1 - b * p_mpa + c * b * p_mpa))
#     return res
# end

function dilation_empirical_function(p_mpa, a, b, c, d)
    res = a * p_mpa /(1 + b*p_mpa) + c * p_mpa + d * p_mpa^2
    return res
end


function dilation_empirical_function_sqerr(X, Y, params)
    sum((dilation_empirical_function.(X, params...) .- Y).^2)
end

function find_dilation_function_params(pressures_mpa::AbstractVector{<:Number}, frac_dilations::AbstractVector{<:Number}, uncertainty_method=nothing)
    start = [1, 1, 1, 1.]
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
    else
        return res
    end
end


