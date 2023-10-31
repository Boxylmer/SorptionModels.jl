struct EmpiricalDilation end

struct EmpiricalDilationModel{T} <: DilationModel
    a::T
    b::T
    c::T
    d::T
end

function EmpiricalDilationModel(a::T, b=zero(T)::T, c=zero(T)::T, d=zero(T)::T) where {T}
    return EmpiricalDilationModel(a, b, c, d)
end

function dilation_empirical_function(p_mpa, a, b=0, c=0, d=0)
    res = @. a * p_mpa + b * p_mpa /(1 + c*p_mpa) + d * p_mpa^2
    return res
end

dilation_empirical_function(p_mpa, m::EmpiricalDilationModel) = dilation_empirical_function(p_mpa, m.a, m.b, m.c, m.d)
predict_dilation(m::EmpiricalDilationModel, pressures_mpa) = dilation_empirical_function(pressures_mpa, m)

"""
    fit_model(::EmpiricalDilation, pressures_mpa, frac_dilations, uncertainty_method=nothing; n_params=3, kwargs...)

Fit the dilation data to an empirical dilation model with no physical meaning. This model is intended to fit any dilation curve given enough parameters. 
The user may specify up to 4 empirical paramteres. Loosely, the use of each parameter is as follows:
- n_params = 1: Linear dilation
- n_params = 2: Naive dual mode like dilation
- n_params = 3: Dual mode like dilation
- n_params = 4: Dual mode like dilation with concave down or tapering behavior permitted

If using `uncertainty_method` (`:JackKnife` or `:Hessian` implemented), each extra parameter will *greatly* increase the resulting uncertainty, as this model is completely empirical.
"""
function fit_model(::EmpiricalDilation, pressures_mpa, frac_dilations, uncertainty_method=nothing; n_params=3, kwargs...)
    start = ones(n_params)
    params = find_dilation_function_params(pressures_mpa, frac_dilations, dilation_empirical_function, uncertainty_method; start, kwargs...)
    return EmpiricalDilationModel(params...)
end
