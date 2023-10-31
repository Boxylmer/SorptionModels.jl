abstract type DilationModel end

"""
    predict_dilation(model::SorptionModel, args...)

Predict fractional dilation with a dilation model. Predictions are generally returned in the format of a vector of values corresponding to an input of a vector of pressures.
"""
function predict_dilation(model, _...)
    t = typeof(model)
    @error "Dilation not implemented for $t"
end

"""
    fit_model(model, pressures_mpa, fractional_dilations, [uncertainty_method])

Fit a dilation model to a set of pressures and fractional dilations. 

Any keyword arguments get passed on to the specific model.

`uncertainty_method` can be any supported method of propagating uncertainty within the model. 
"""
function fit_model(model::DilationModel, _...)
    t = typeof(model)
    @error "Fitting not implemented for $t"
end


dilation_function_sqerr(X, Y, params, func::Function) = sum((func.(X, params...) .- Y).^2)

function find_dilation_function_params(pressures_mpa, frac_dilations, func::Function, uncertainty_method=nothing; start=ones(n_params))
    ps, fs = strip_measurement_to_value(pressures_mpa), strip_measurement_to_value(frac_dilations) 
    obj = x -> dilation_function_sqerr(ps, fs, x, func)
    n = length(pressures_mpa)
    res = Optim.optimize(obj, start, BFGS()).minimizer

    if uncertainty_method == :Hessian
        model_uncertainty = rss_minimizer_standard_errors(obj, res, n)
        return res .± model_uncertainty

    elseif uncertainty_method == :JackKnife
        data = collect(zip(pressures_mpa, frac_dilations))
        fit_function(data) = find_dilation_function_params(collect.(collect(zip(data...)))..., func; start = res)
        model_uncertainty = jackknife_uncertainty(fit_function, data)
        return res .± model_uncertainty

    elseif isnothing(uncertainty_method)
        return res
    else
        throw(ArgumentError("$uncertainty_method not supported. See available error propagation methods in DilationModels."))
    end
end

# todo docs
"Get the derivative of a dilation vs pressure curve at a given pressure."
predict_dilation_derivative(dilation_model::DilationModel, pressures_mpa) = ForwardDiff.derivative.(x -> predict_dilation(dilation_model, x), pressures_mpa)