struct FloryHugginsDualMode end

struct FloryHugginsDualModeModel{CHIT, BT, KDT} <: SorptionModel
    ch::CHIT  # units: CC/CC
    b::BT    # units: None
    chi::KDT  # units: None
    poldens::PDT # units: g/cm3
end

function Base.show(io::IO, obj::FloryHugginsDualModeModel)
    ch = obj.ch; b = obj.b; kd = obj.chi
    print(io, "Flory Huggins Dual Mode: CH' = $ch, b = $b, χ = $kd")
end

"""
    infinite_dilution_solubility(dm::FloryHugginsDualModeModel)
Get the infinite dilution solubility coefficient in (**(CC/CC)**) (activity is unitless)
"""
infinite_dilution_solubility(dm::FloryHugginsDualModeModel) = nothing # TODO


"""
    a_predict_concentration(::FloryHugginsDualModeModel, activity::Number)
Predict the concentration of a penetrant given a flory huggins dual mode model and activity.
"""
function a_predict_concentration(fhdm::FloryHugginsDualModeModel, activity::Number)
    fh_mode = nothing # TODO
    return fh_mode + dm.ch * dm.b * activity / (1 + fhdm.b * activity)
end

function MembraneBase.rss(dm::FloryHugginsDualMode, isotherm::IsothermData)
    if isotherm.num_components == 1
        if dm.use_fugacity
            predictions = predict_concentration(dm, fugacities(isotherm; component=1))  
        else
            predictions = predict_concentration(dm, partial_pressures(isotherm, component=1))
        end
        return MembraneBase.rss(concentration(isotherm, component=1), predictions)
    end

    # no implementation yet for multiple components
    throw(DimensionMismatch("The isotherm given has more than one component, this function only works for pure isotherms")) 
end

# function MembraneBase.rss(dm::DualModeModel, pres_or_fugs::AbstractVector, concentrations::AbstractVector)
#     return rss(predict_concentration(dm, pres_or_fugs), concentrations)
# end


# function fit_dualmode_model(isotherm::IsothermData; uncertainty_method=nothing, use_fugacity=false, apply_weights=false)
#     # see if isotherm is only a single component
#     if isotherm.num_components != 1
#         throw(ErrorException("The isotherm given has more than one component, this function only works for pure isotherms"))
#     end

#     if !apply_weights
#         used_isotherm = strip_measurement_to_value(isotherm)
#     else
#         used_isotherm = isotherm
#     end

#     target = function(ch_b_kd)
#         dmm = DualModeModel(ch_b_kd...; use_fugacity)
#         err = rss(dmm, used_isotherm)
#         if typeof(err) <: Measurement err = err.val end  # handle measurement types (we don't need them where we're going!)
#         return err
#     end
#     lower = [0., 0., 0.]
#     upper = [Inf, Inf, Inf]
#     res = Optim.optimize(target, lower, upper, [0.5, 0.5, 0.5], Fminbox(BFGS()))
    
#     if use_fugacity
#         pressure_function = fugacities
#         resampled_fitting_function = DualModeHelperFunctions.resampled_set_fitting_wrapper_fugacity
#     else 
#         pressure_function = partial_pressures
#         resampled_fitting_function = DualModeHelperFunctions.resampled_set_fitting_wrapper_pressure
#     end
    
#     if !isnothing(uncertainty_method)
#         ps = pressure_function(used_isotherm; component=1)
#         cs = concentration(used_isotherm; component=1)
#         data = collect(zip(ps, cs))
#     end

#     if uncertainty_method == :JackKnife
#         corresponding_uncertainties = jackknife_uncertainty(resampled_fitting_function, data)
#         uncertain_parameters = [Optim.minimizer(res)[i] ± corresponding_uncertainties[i] for i in 1:length(corresponding_uncertainties)]
#         optimized_model = DualModeModel(uncertain_parameters...; use_fugacity)
#     elseif uncertainty_method == :Bootstrap
#         corresponding_uncertainties = bootstrap_uncertainty(resampled_fitting_function, data)  
#         uncertain_parameters = [Optim.minimizer(res)[i] ± corresponding_uncertainties[i] for i in 1:length(corresponding_uncertainties)]
#         optimized_model = DualModeModel(uncertain_parameters...; use_fugacity)  
#     elseif isnothing(uncertainty_method)
#         optimized_model = DualModeModel(Optim.minimizer(res)...; use_fugacity)
#     else
#         throw(ArgumentError("Invalid uncertainty_method: " * string(uncertainty_method)))
#     end
#     return optimized_model
# end

# """
#     fit_model(DualMode(), isotherm::IsothermData, [uncertainty_method=nothing], [apply_weights=false], [use_fugacity=false])
# Fit the dual mode model to the pressures and concentrations present in the isotherm. 

# Options
# - For determining the uncertainty of the model parameters, the `:JackKnife`, and `:Bootstrap` methods are available. 
# - `apply_weights` will use a weighted nonlinear regression method to solve the parameters, given that `Measurement` types are used somewhere in the data. 
# - `use_fugacity` will fit the model to fugacities instead of pressures (they should be present in the isotherm data).  
# """
# fit_model(::DualMode, isotherm::IsothermData; kwargs...) = fit_dualmode_model(isotherm; kwargs...)

# function MembraneBase.strip_measurement_to_value(model::DualModeModel)
#     return DualModeModel(
#         strip_measurement_to_value(model.ch),
#         strip_measurement_to_value(model.b),
#         strip_measurement_to_value(model.kd),
#         model.use_fugacity
#         )
# end