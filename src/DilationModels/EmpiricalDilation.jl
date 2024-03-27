struct EmpiricalDilation{N} end
EmpiricalDilation(n=3) = EmpiricalDilation{n}()
n_adjustable_params(ed::EmpiricalDilation) = typeof(ed).parameters[1]

struct EmpiricalDilationModel{T} <: DilationModel
    a::T
    b::T
    c::T
    d::T
    e::T
end

function EmpiricalDilationModel(a::T, b=zero(T)::T, c=zero(T)::T, d=zero(T)::T, e=zero(T)::T) where {T}
    return EmpiricalDilationModel(a, b, c, d, e)
end

function dilation_empirical_function(p_mpa, a, b=0, c=0, d=0, e=0)
    res = @. a * p_mpa + b * p_mpa /(1 + c*p_mpa) + d * exp(e * p_mpa)
    return res
end

dilation_empirical_function(p_mpa, m::EmpiricalDilationModel) = dilation_empirical_function(p_mpa, m.a, m.b, m.c, m.d, m.e)
predict_dilation(m::EmpiricalDilationModel, pressures_mpa) = dilation_empirical_function(pressures_mpa, m)

"""
    fit_model(::EmpiricalDilation, pressures_mpa, frac_dilations, uncertainty_method=nothing; n_params=3, kwargs...)

Fit the dilation data to an empirical dilation model with no physical meaning. This model is intended to fit any dilation curve given enough parameters. 
The user may specify up to 4 empirical paramteres. Loosely, the use of each parameter is as follows:
- n_params = 1: Linear dilation
- n_params = 2-3: Unscaled(2) dual mode like (3) dilation (3=default)
- n_params = 4-5: Dual mode like dilation with concave up (4) or tapering behavior permitted with scaling (5)

If using `uncertainty_method` (`:JackKnife` or `:Hessian` implemented), each extra parameter will *greatly* increase the resulting uncertainty, as this model is completely empirical.
"""
function fit_model(ed::EmpiricalDilation, pressures_mpa, frac_dilations; kwargs...)
    start = ones(n_adjustable_params(ed))
    params = find_dilation_function_params(pressures_mpa, frac_dilations, dilation_empirical_function, start; kwargs...)
    return EmpiricalDilationModel(params...)
end
