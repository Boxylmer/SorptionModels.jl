struct GAB end

struct GABModel{CPT, KT, AT, PCFT, ACFT} <: SorptionModel
    cp::CPT  # units: CC/CC
    k::KT    # units: None
    a::AT    # units:   None 
    pressure_conversion_function::PCFT  # MPa to activity
    activity_conversion_function::ACFT  # Activity to MPa
end

GABModel(cp, k, a) = GABModel(cp, k, a, missing, missing)
    

module GABHelperFunctions
    using SorptionModels
    using MembraneBase
    function resampled_set_fitting_wrapper(resampled_dataset)
        # here dataset is in the format of [(a1, c1), (a2, c2), ...] This entire thing is quick and *dirty*
        
        # return resampled dataset to [a1, a2, ...], [c1, c2, ...] format
        _formatted_dataset = collect(zip(resampled_dataset...))
        activity_vector = collect(_formatted_dataset[1])
        concentration_vector = collect(_formatted_dataset[2])

        # acquire and vectorize fitting
        gmfitting = fit_gab_model(activity_vector, concentration_vector)  # collect converts iterables to arrays
        return [gmfitting.cp, gmfitting.k, gmfitting.a]
    end

end #  end helper functions

"""
    predict_concentration(gm::GABModel, pressure::Number)
Predict a concentration given a pressure. Note that GAB is fit to activity, so you will need to specify a `pressure_conversion_function(pressure::Number) -> activity::Number` when fitting the model.
See `fit_gab_model`
"""
predict_concentration(gm::GABModel, pressure::Number) = a_predict_concentration(gm, gm.pressure_conversion_function(pressure))
predict_concentration(gm::GABModel, pressures::AbstractVector) = a_predict_concentration(gm, gm.pressure_conversion_function.(pressures))


"""
    a_predict_concentration(gm::GABModel, activity::Number)
Predict concentration based on a GAB model and given activity.
Mixed GAB models aren't implemented due to the empirical nature of the model (and it's lack of predictivity)

"""
function a_predict_concentration(gm::GABModel, activity::Number)
    (gm.cp * gm.k * gm.a * activity) / ((1 - gm.k*activity) * (1 - gm.k*activity + gm.k*gm.a*activity))
end

"""
    a_predict_concentration(gm::GABModel, activities::AbstractVector)
Returns a vector of concentrations predicted by the GAB model for a corresponding vector of activities.
"""
function a_predict_concentration(gm::GABModel, activities::AbstractVector)
    predictions = [predict_concentration(gm, activity) for activity in activities]
    return predictions
end

function MembraneBase.rss(gm::GABModel, activities::AbstractVector, concentrations::AbstractVector)
    if length(activities) == length(concentrations)
        predictions = a_predict_concentration(gm, activities)
        return MembraneBase.rss(concentrations, predictions)        
    end

    # no implementation for multiple components
    throw(ErrorException("The isotherm given has more than one component, this function only works for pure isotherms")) 
end

"""
    fit_gab_model(activities::AbstractVector, concentrations::AbstractVector; uncertainty_method=nothing)
Fit a set of activities and corresponding concentrations (**CC/CC**) to the GAB model. 

For determining the uncertainty of the model parameters, the `:JackKnife`, and `:Bootstrap` methods are available. 

- `pressure_conversion_function(pressure::Number) -> activity::Number`: if set, will allow pressures to automatically be converted to activities when needed. (MPa) 
- `activity_conversion_function(activity::Number) -> pressure::Number`: if set, will allow activities to automatically be converted to pressures when needed. (MPa)
"""
function fit_gab_model(activities::AbstractVector, concentrations::AbstractVector; 
    uncertainty_method=nothing, apply_weights=false, pressure_conversion_function=missing, activity_conversion_function=missing)

    if length(activities) != length(concentrations)
        throw(DimensionMismatch("The concentrations and activities given don't match in length, this function only works for equivalent length vectors of numbers"))
    end

    if !apply_weights
        applied_activities = strip_measurement_to_value(activities)
        applied_concs = strip_measurement_to_value(concentrations)
    else
        applied_activities = activities
        applied_concs = concentrations
    end

    # target = function(cp_k_a)
    #     gm = GABModel(cp_k_a...)
    #     err = rss(gm, applied_activities, applied_concs)
    #     # if typeof(err) <: Measurement err = err.val end  # handle measurement types (we don't need them where we're going!)
    #     return err
    # end
    target = _make_gab_target(applied_activities, applied_concs)

    # res = Optim.optimize(target, [1., 1., 1.]; autodiff = :forward)
    # lower = [0., 0., 0.]
    # upper = [Inf, Inf, Inf]
    # res = Optim.optimize(target, lower, upper, [1., 1., 1. ], Fminbox(BFGS()); autodiff = :forward)
    init_params = _get_initial_conditioned_gab_params(applied_activities, applied_concs)
    res = Optim.optimize(target, init_params, Newton(); autodiff = :forward,)


    # fit_params = _uncondition_gab_guess(Optim.minimizer(res))
    fit_params = Optim.minimizer(res)

    # todo this can be made much more readable by moving corresponding_uncertainties and uncertain_parameters out of the ifs and introducing dummy values in the else case.
    if uncertainty_method == :JackKnife
        data = collect(zip(applied_activities, applied_concs))
        corresponding_uncertainties = jackknife_uncertainty(GABHelperFunctions.resampled_set_fitting_wrapper, data)
        final_params = [fit_params[i] ± corresponding_uncertainties[i] for i in 1:length(corresponding_uncertainties)]
    elseif uncertainty_method == :Bootstrap
        data = collect(zip(applied_activities, applied_concs))  # Isotherm.dataset(isotherm)
        corresponding_uncertainties = bootstrap_uncertainty(GABHelperFunctions.resampled_set_fitting_wrapper, data)  
        final_params = [fit_params[i] ± corresponding_uncertainties[i] for i in 1:length(corresponding_uncertainties)]
    elseif isnothing(uncertainty_method)
        final_params = fit_params
    else 
        throw(ArgumentError("Invalid uncertainty_method: " * string(uncertainty_method)))
    end
    optimized_model = GABModel(final_params..., pressure_conversion_function, activity_conversion_function)
    return optimized_model
end

# function _condition_gab_guess(guess)
#     return [guess[1], log(guess[2]), log(guess[3])]
# end

# _uncondition_gab_guess(guess) = 10 .^ guess

function _get_initial_conditioned_gab_params(applied_activities, applied_concs)
    cp = maximum(applied_concs)
    target = _make_gab_presolved_target(applied_activities, applied_concs, cp)
    res = Optim.optimize(target, [0.01, 100], Newton(); autodiff = :forward)
    cp_k_a = [cp, Optim.minimizer(res)...]
    return cp_k_a
end

function _make_gab_presolved_target(applied_activities, applied_concs, cp_val)
    target = function(k_a)
        # k_a = _uncondition_gab_guess(k_a_conditioned_guess)
        cp_k_a = [cp_val, k_a...]
        gm = GABModel(cp_k_a...)
        err = rss(gm, applied_activities, applied_concs)
        # if typeof(err) <: Measurement err = err.val end  # handle measurement types (we don't need them where we're going!)
        return err
    end
    return target
end

function _make_gab_target(applied_activities, applied_concs)
    target = function(cp_k_a)
        # cp_k_a = _uncondition_gab_guess(cp_k_a_conditioned_guess)
        # @show cp_k_a
        gm = GABModel(cp_k_a...)
        err = rss(gm, applied_activities, applied_concs)
        # if typeof(err) <: Measurement err = err.val end  # handle measurement types (we don't need them where we're going!)
        return err
    end
    return target
end


"""
    fit_gab_model(isotherm::IsothermData; kwargs...)
Fit the GAB model to the concentrations and activities present in an isotherm. 
Right now, it is assumed that the isotherm has only one component, so the model fits only to the first component. 
# todo: if a multicomponent isotherm is provided, return a vector of GAB fittings. 


"""
function fit_gab_model(isotherm::IsothermData; kwargs...)
    acts = activities(isotherm; component=1)
    concs = concentration(isotherm; component=1)
    if isnothing(acts)
        throw(MissingException("Isotherm has no activities."))
    elseif isnothing(concs)
        throw(MissingException("Isotherm has no concentrations.")) 
    end

    return fit_gab_model(acts, concs; kwargs...)
end

fit_model(::GAB, isotherm::IsothermData; kwargs...) = fit_gab_model(isotherm; kwargs...)
fit_model(::GAB, activities::AbstractVector, concentrations::AbstractVector; kwargs...) = fit_gab_model(activities, concentrations; kwargs...)