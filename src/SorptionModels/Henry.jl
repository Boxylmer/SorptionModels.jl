struct Henry end

struct HenryModel{KDT} <: SorptionModel
    kd::KDT  # units: CC/CC / MPa
    use_fugacity::Bool
end

function HenryModel(kd; use_fugacity=false)
    return HenryModel(kd, use_fugacity)
end

function Base.show(io::IO, obj::HenryModel)
    kd = obj.kd
    using_fug = obj.use_fugacity ? "using fugacity." : "using pressures."
    print(io, "Henry's Law: kd = $kd $using_fug")
end

"""
    infinite_dilution_solubility(dm::HenryModel)
Get the infinite dilution solubility coefficient in (**(CC/CC) / MPa**))
"""
infinite_dilution_solubility(dm::Henry) = dm.kd


"""
    predict_concentration(::HenryModel, pressure_mpa::Number)
Predict the concentration of a penetrant given a henry model and pressure (**MPa**).
"""
function predict_concentration(hm::HenryModel, pressure_mpa::Number)
    hm.kd * pressure_mpa
end



"""
Predict a pressure (**MPa**) given a concentration (**CC(STP)/CC(Polymer)**) according to Henry's Law.
"""
function predict_pressure(dm::HenryModel, concentrations_cc_cc::Number)
    k = dm.kd
    c = concentrations_cc_cc
    return c / k
end

function MembraneBase.rss(hm::HenryModel, isotherm::IsothermData)
    if isotherm.num_components == 1
        if hm.use_fugacity
            predictions = predict_concentration.(hm, fugacities(isotherm; component=1))  
        else
            predictions = predict_concentration.(hm, partial_pressures(isotherm, component=1))
        end
        return MembraneBase.rss(concentration(isotherm, component=1), predictions)
    end

    # no implementation yet for multiple components
    throw(DimensionMismatch("The isotherm given has more than one component, this function only works for pure isotherms")) 
end

function MembraneBase.rss(hm::HenryModel, pres_or_fugs::AbstractVector, concentrations::AbstractVector)
    return rss(predict_concentration.(hm, pres_or_fugs), concentrations) # TODO can refactor to a general sorptionmodel
end

# TODO This needs to just be linest with intercept = 0, but right now I'm in a hurry to get permeation things working for helium data
function fit_henry_model(isotherm::IsothermData; uncertainty_method=nothing, use_fugacity=false, apply_weights=false)
    # see if isotherm is only a single component
    if isotherm.num_components != 1
        throw(ErrorException("The isotherm given has more than one component, this function only works for pure isotherms"))
    end

    if !apply_weights
        used_isotherm = strip_measurement_to_value(isotherm)
    else
        used_isotherm = isotherm
    end

    target = function(kd)
        dmm = HenryModel(kd...; use_fugacity)
        err = rss(dmm, used_isotherm)
        if typeof(err) <: Measurement err = err.val end  
        return err
    end
    lower = [0., ]
    upper = [Inf]
    res = Optim.optimize(target, lower, upper, [0.5,], Fminbox(BFGS())).minimizer
    

    if !isnothing(uncertainty_method)
        throw(ArgumentError("This model supports no uncertainty methods"))
    end
    return HenryModel(res..., use_fugacity)
end

"""
    fit_model(Henry(), isotherm::IsothermData, [uncertainty_method=nothing], [apply_weights=false], [use_fugacity=false])
Fit Henry's Law to the pressures and concentrations present in the isotherm. 

"""
fit_model(::Henry, isotherm::IsothermData; kwargs...) = fit_henry_model(isotherm; kwargs...)


function MembraneBase.strip_measurement_to_value(model::HenryModel)
    return HenryModel(
        strip_measurement_to_value(model.kd),
        model.use_fugacity
    )
end

