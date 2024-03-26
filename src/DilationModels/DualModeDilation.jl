struct DualModeDilation end
n_adjustable_params(::DualModeDilation) = 2

struct DualModeDilationModel{T} <: DilationModel
    vd::T  # cm3(pol) / cm3(stp)
    f::T   # unitless
    ch::T  # cc/cc
    b::T   # 1/MPa
    kd::T  # cc/cc / MPa
end

function DualModeDilationModel(vd, f, m::DualModeModel)
    return DualModeDilationModel(promote(vd, f, m.ch, m.b, m.kd)...)
end

dualmode_dilation_function(p_mpa, vd, f, ch, b, kd) = @. vd * (kd * p_mpa + f*ch*b*p_mpa / (1 + b*p_mpa))
dualmode_dilation_function(p_mpa, vd, f, m::DualModeModel) = dualmode_dilation_function(p_mpa, vd, f, m.ch, m.b, m.kd)
dualmode_dilation_function(p_mpa, m::DualModeDilationModel) = dualmode_dilation_function(p_mpa, m.vd, m.f, m.ch, m.b, m.kd)
predict_dilation(m::DualModeDilationModel, pressures_mpa) = dualmode_dilation_function(pressures_mpa, m)

"""
    fit_model(::DualModeDilation, pressures_mpa, frac_dilations, uncertainty_method=nothing; kwargs...)

Fit the dilation data to an dual mode dilation model, using a dual mode model previoulsy fit to concentration and pressure data. 


Dilation is described by two additional parameters, ``V_D`` and ``f`` using the following model:

``{\\Delta}{V} / V_0 = V_D (k_D \\cdot p + \\frac{f C_{H}^{\'}bp}{1+bp})``
- ``V_D`` represents the effective condensed penetrant molar volume.
- ``f`` represents an empirical fit describing how much the langmuir mode participates in dilation.


`J.D. Moon, M. Galizia, H. Borjigin, R. Liu, J.S. Riffle, B.D. Freeman, D.R. Paul, Water Vapor Sorption, Diffusion, and Dilation in Polybenzimidazoles, Macromolecules. 51 (2018) 7197–7208. https://doi.org/10.1021/acs.macromol.8b01659.`
"""
function fit_model(::DualModeDilation, pressures_mpa, frac_dilations, model::DualModeModel; kwargs...)
    if model.use_fugacity
        throw(ArgumentError("Dual mode model must be fit to pressure, not fugacity, activity, or other potentials."))
    end
    start = ones(n_adjustable_params(DualModeDilation()))
    pure_model = strip_measurement_to_value(model)
    func(p, vd, f) = dualmode_dilation_function(p, vd, f, pure_model)
    params = find_dilation_function_params(pressures_mpa, frac_dilations, func, start; kwargs...)
    return DualModeDilationModel(params..., model)
end
