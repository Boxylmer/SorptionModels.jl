struct GAB end

struct GABModel{CPT, KT, AT} <: SorptionModel
    cp::CPT  # units: CC/CC
    k::KT    # units: None
    a::AT    # units:   None 
end

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
    predict_concentration(gm::GABModel, activity::Number)
Predict concentration based on a GAB model and given activity.
Mixed GAB models aren't implemented due to the empirical nature of the model (and it's lack of predictivity)

"""
function predict_concentration(gm::GABModel, activity::Number)
    (gm.cp * gm.k * gm.a * activity) / ((1 - gm.k*activity) * (1 - gm.k*activity + gm.k*gm.a*activity))
end

"""
    predict_concentration(gm::GABModel, activities::AbstractVector)
Returns a vector of concentrations predicted by the GAB model for a corresponding vector of activities.
"""
function predict_concentration(gm::GABModel, activities::AbstractVector)
    predictions = [predict_concentration(gm, activity) for activity in activities]
    return predictions
end

function rss(gm::GABModel, activities::AbstractVector, concentrations::AbstractVector)
    if length(activities) == length(concentrations)
        predictions = predict_concentration(gm, activities)
        return MembraneBase.rss(concentrations, predictions)        
    end

    # no implementation for multiple components
    throw(ErrorException("The isotherm given has more than one component, this function only works for pure isotherms")) 
end

"""
    fit_gab_model(activities::AbstractVector, concentrations::AbstractVector; uncertainty_method=nothing)
Fit a set of activities and corresponding concentrations (**CC/CC**) to the GAB model. 

For determining the uncertainty of the model parameters, the `:JackKnife`, and `:Bootstrap` methods are available. 
"""
function fit_gab_model(activities::AbstractVector, concentrations::AbstractVector; uncertainty_method=nothing)
    if length(activities) != length(concentrations)
        throw(DimensionMismatch("The concentrations and activities given don't match in length, this function only works for equivalent length vectors of numbers"))
    end
    target = function(cp_k_a)
        gm = GABModel(cp_k_a...)
        err = rss(gm, activities, concentrations)
        if typeof(err) <: Measurement err = err.val end  # handle measurement types (we don't need them where we're going!)
        return err
    end
    
    res = Optim.optimize(target, [1., 1., 1.]; autodiff = :forward)
    
    # todo this can be made much more readable by moving corresponding_uncertainties and uncertain_parameters out of the ifs and introducing dummy values in the else case.
    if uncertainty_method == :JackKnife
        data = collect(zip(activities, concentrations))
        corresponding_uncertainties = jackknife_uncertainty(GABHelperFunctions.resampled_set_fitting_wrapper, data)
        uncertain_parameters = [Optim.minimizer(res)[i] ± corresponding_uncertainties[i] for i in 1:length(corresponding_uncertainties)]
        optimized_model = GABModel(uncertain_parameters...)
    elseif uncertainty_method == :Bootstrap
        data = collect(zip(activities, concentrations))  # Isotherm.dataset(isotherm)
        corresponding_uncertainties = bootstrap_uncertainty(GABHelperFunctions.resampled_set_fitting_wrapper, data)  
        uncertain_parameters = [Optim.minimizer(res)[i] ± corresponding_uncertainties[i] for i in 1:length(corresponding_uncertainties)]
        optimized_model = GABModel(uncertain_parameters...)  
    elseif isnothing(uncertainty_method)
        optimized_model = GABModel(Optim.minimizer(res)...)
    else 
        throw(ArgumentError("Invalid uncertainty_method: " * string(uncertainty_method)))
    end
    return optimized_model
end

"""
    fit_gab_model(isotherm::IsothermData; kwargs...)
Fit the GAB model to the concentrations and activities present in an isotherm. 
Right now, it is assumed that the isotherm has only one component, so the model fits only to the first component. 
# todo: if a multicomponent isotherm is provided, return a vector of GAB fittings. 
"""
fit_gab_model(isotherm::IsothermData; kwargs...) = fit_gab_model(activities(isotherm; component=1), concentration(isotherm; component=1); kwargs...)

fit_model(::GAB, isotherm::IsothermData; kwargs...) = fit_gab_model(isotherm; kwargs...)
