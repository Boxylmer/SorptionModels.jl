struct DualModeDilation end

struct DualModeDilationModel{T} <: DilationModel
    vd::T
    f::T
    ch::T
    b::T
    kd::T
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
# todo include model math here
"""
function fit_model(::DualModeDilation, pressures_mpa, frac_dilations, model::DualModeModel, uncertainty_method=nothing; kwargs...)
    if model.use_fugacity
        throw(ArgumentError("Dual mode model must be fit to pressure, not fugacity, activity, or other potentials."))
    end
    start = ones(2)
    pure_model = strip_measurement_to_value(model)
    func(p, vd, f) = dualmode_dilation_function(p, vd, f, pure_model)
    params = find_dilation_function_params(pressures_mpa, frac_dilations, func, uncertainty_method; start, kwargs...)
    return DualModeDilationModel(params..., model)
end
