abstract type TransientSorptionModel end 

"""
    predict_sorption(sorptionmodel::TransientSorptionModel, time_seconds::AbstractVector{<:Number})
Predict the amount of dimensionless sorption that has occurred at time `time_seconds` according to the transient sorption model.

!!! note
    You can also specify a single time instead of a vector of times. 
    When predicting this way, you can additionally use the keyword argument `iters`, which will change the default number of iterations used to solve the transient sorption prediction. 
    It is not recommended to modify this default value. 
"""
function predict_sorption(sorptionmodel::TransientSorptionModel, time_seconds::AbstractVector{<:Number})
    return [predict_sorption(sorptionmodel, t) for t in time_seconds]
end

"""
    rss(sorptionmodel::TransientSorptionModel, step_data::TransientStepData)
Determine how well a model predicts transient sorption by the sum of squared residuals.
"""
function rss(sorptionmodel::TransientSorptionModel, step_data::TransientStepData)
    sorption_pred = predict_sorption(sorptionmodel, step_data.time)
    sorption_err = sum((sorption_pred - step_data.dimensionlesssorption).^2)
    return sorption_err
end

"""
    fit_transient_sorption_model(
        step_data::TransientStepData, model_symbol; custom_initial_params=nothing, interpolation_method=nothing, interpolation_datapoints=1000, 
        uncertainty_method=nothing, num_uncertainty_resamples=20, resampling_mode=false)

# Arguments
- `step_data::TransientStepData`: See [`TransientStepData`](@ref)
- `model_symbol`: Type of transient sorption model to fit the `step_data` to. Available options are:
    - `:FickianSorptionModel`
    - `:BerensHopfenbergSorptionModel`
    - `:ModifiedBerensHopfenbergSorptionModel`
# Optional Arguments
- `custom_initial_params::Vector`: Vector of initial guesses to start with. Must match the number of parameters in the model.
    - These parameters are *linearized* in their parameter space. Make sure you understand the linearization of the model parameters before you try to specify custom initial guesses. 
- `interpolation_method`: See `resample` in [`TransientStepData`](@ref).
- `interpolation_datapoints::Integer`: Number of data points to interpolate to, see above. 
- `uncertainty_method`: currently only supports `nothing` and the `:Bootstrap` method of determining uncertainty.
- `num_uncertainty_resamples`: Number of times to resample for determining the uncertainty of the model fit (if applicable).
    - e.g., for `:Bootstrap`, this would be the number of times the data is randomly resampled and refit.
- `resampling_mode`: For internal use with Optim.jl and resampling. Do not modify this value. 
    - This will return the linearized parameter fittings in lieu of the normal model objects. 
"""
function fit_transient_sorption_model(
    step_data::TransientStepData, model_symbol; custom_initial_params=nothing, interpolation_method=nothing, interpolation_datapoints=1000, 
    uncertainty_method=nothing, num_uncertainty_resamples=20, resampling_mode=false)
    
    # resample data if necessary
    if !isnothing(interpolation_method)
        _step_data = strip_measurement_to_value(resample(step_data, interpolation_datapoints, interpolation_method))
    else 
        _step_data = strip_measurement_to_value(step_data) 
    end
    # prepare fitting model and parameters
    if model_symbol == :FickianSorptionModel
        ub = [1.1, log(1.)]
        lb = [0., -50.]
        initial_params = [1., log(0.001)]
        model = FickianSorptionModel
    elseif model_symbol == :BerensHopfenbergSorptionModel
        ub = [1.1, log(1.), 1., log(1.)]
        lb = [0., -50., 0., -50.]
        # initial_params = [0.9, -3, 0.1, -10]
        model = BerensHopfenbergSorptionModel
        if isnothing(custom_initial_params)
            initial_fickan_fit = linearize_model(fit_transient_sorption_model(_step_data, :FickianSorptionModel))
            initial_params = [initial_fickan_fit.m_f, initial_fickan_fit.k_f, 0.1, log(0.0001)]
        end
    elseif model_symbol == :ModifiedBerensHopfenbergSorptionModel
        ub = [1.1, log(1.), 1., log(1.), log(10)]
        lb = [0., -50., 0., -50., -30.]
        # initial_params = [0.9, -3, 0.1, -10, -5]
        model = ModifiedBerensHopfenbergSorptionModel
        if isnothing(custom_initial_params)
            initial_bh_fit = linearize_model(fit_transient_sorption_model(_step_data, :BerensHopfenbergSorptionModel))
            initial_params = [initial_bh_fit.m_f, initial_bh_fit.k_f, initial_bh_fit.m_r, initial_bh_fit.k_r, 1.]
        end
    else
        throw(ArgumentError("Specified model symbol not found, try \":BerensHopfenbergSorptionModel\"?"))   
    end

    if !isnothing(custom_initial_params)  # use given initial parameters if specified
        initial_params = custom_initial_params
    end

    # actually do the fitting
    target_function = function(params)
        linearized_model = model(params..., true)
        unlinearized_model = unlinearize_model(linearized_model)
        return rss(unlinearized_model, _step_data)
    end
    if resampling_mode
        iterations= 10
        time_limit = 1
    else
        iterations= 100
        time_limit = 120
    end

    # res = Optim.optimize(target_function, lb, ub, initial_params, Fminbox(BFGS()), 
    #     Optim.Options(allow_f_increases=true, iterations=iterations, time_limit=time_limit); autodiff=:forward)
    
    res = Optim.optimize(target_function, lb, ub, initial_params, Fminbox(BFGS()), 
        Optim.Options(allow_f_increases=true, iterations=iterations, time_limit=time_limit); autodiff=:forward)

    fitted_params = Optim.minimizer(res)
    # deal with fittings
    "Specify a dataset [(time, sorption)...] and get back a vector of parameters."
    function fitting_uncertainty_wrapper(full_dataset)
        transient_step = TransientStepData(full_dataset)
        return fit_transient_sorption_model(transient_step, model_symbol; custom_initial_params=fitted_params, resampling_mode=true, uncertainty_method=nothing)
    end

    if uncertainty_method == :Bootstrap
        errors = bootstrap_uncertainty(
            fitting_uncertainty_wrapper, dataset(_step_data); nsamples=num_uncertainty_resamples
        )
    
    elseif uncertainty_method == :JackKnife
        throw(ErrorException("Jackknife isn't working for this set of models just yet"))
        errors = jackknife_uncertainty(fitting_uncertainty_wrapper, dataset(_step_data))

    elseif uncertainty_method == :Hessian
        throw(ErrorException("Inverse Hessian error propagation isn't working for this set of models just yet"))
    elseif isnothing(uncertainty_method)
        errors = nothing
    else
        throw(ArgumentError("Uncertainty method should be a symbol, e.g., :JackKnife or :Bootstrap"))
    end

    # if we did uncertainty calculations, lets make measurement types instead
    if !isnothing(errors)
        fitted_params = fitted_params .± errors
    end 

    # if we're doing an uncertainty_method like :Bootstrap or :JackKnife, 
    #   we need to return the vector of results instead of a finished model
    if resampling_mode
        return fitted_params
    else
        return unlinearize_model(model(fitted_params..., true))
    end
end

"""
    get_diffusivity(model::TransientSorptionModel, semi_thickness_cm::Number)
Get the diffusivity of the polymer described by the transient sorption model. 
- The semi-thickness (half of the polymer thickness, in **cm**) of the polymer used to collect the sorption data is required.
"""
function get_diffusivity(model::TransientSorptionModel, semi_thickness_cm::Number)
    return model.k_f * (2 * semi_thickness_cm)^2 / pi^2
end
get_diffusivity(model::TransientSorptionModel, semi_thickness_cm::Missing) = missing  # isn't this redundant?