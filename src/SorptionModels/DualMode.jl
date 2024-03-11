struct DualMode end

struct DualModeModel{CHT, BT, KDT} <: SorptionModel
    ch::CHT  # units: CC/CC
    b::BT    # units: 1 / MPa
    kd::KDT  # units: CC/CC / MPa
    use_fugacity::Bool
end

function DualModeModel(ch, b, kd; use_fugacity=false)
    return DualModeModel(ch, b, kd, use_fugacity)
end

function Base.show(io::IO, obj::DualModeModel)
    ch = obj.ch; b = obj.b; kd = obj.kd
    using_fug = obj.use_fugacity ? "using fugacity." : "using pressures."
    print(io, "Dual Mode: CH' = $ch, b = $b, kd = $kd $using_fug")
end

module DualModeHelperFunctions
    using SorptionModels
    using MembraneBase
    function resampled_set_fitting_wrapper_pressure(dataset)
        pps, concs = collect(zip(dataset...))
        isotherm = IsothermData(; partial_pressures_mpa = pps, concentrations_cc = concs)
        dmfitting = fit_dualmode_model(isotherm)
        return [dmfitting.ch, dmfitting.b, dmfitting.kd]
    end

    function resampled_set_fitting_wrapper_fugacity(dataset)
        fugs, concs = collect(zip(dataset...))
        isotherm = IsothermData(; fugacities_mpa = fugs, concentrations_cc = concs)
        dmfitting = fit_dualmode_model(isotherm; use_fugacity = true)
        return [dmfitting.ch, dmfitting.b, dmfitting.kd]
    end

end

function langmuir_mode_concentration(dm::DualModeModel, pressure_mpa::Number)
    dm.ch * dm.b * pressure_mpa / (1 + dm.b * pressure_mpa)
end
langmuir_mode_concentration(dm::DualModeModel, pressures_mpa::AbstractVector{<:Number}) = [langmuir_mode_concentration(dm, p) for p in pressures_mpa]

function henry_mode_concentration(dm::DualModeModel, pressure_mpa::Number)
    dm.kd * pressure_mpa
end 
henry_mode_concentration(dm::DualModeModel, pressures_mpa::AbstractVector{<:Number}) = [henry_mode_concentration(dm, p) for p in pressures_mpa]

"""
    infinite_dilution_solubility(dm::DualModeModel)
Get the infinite dilution solubility coefficient in (**(CC/CC) / MPa**))
"""
infinite_dilution_solubility(dm::DualModeModel) = dm.kd + dm.ch*dm.b


"""
    predict_concentration(::DualModeModel, pressure_mpa::Number)
Predict the concentration of a penetrant given a dual mode model and pressure (**MPa**).
"""
function predict_concentration(dm::DualModeModel, pressure_mpa::Number)
    henry_mode_concentration(dm, pressure_mpa) + langmuir_mode_concentration(dm, pressure_mpa)
end

"""
    predict_concentration(::DualModeModel, pressures_mpa::AbstractVector)
Predict a vector of concentrations for a pure dual mode model given a corresponding vector of pressures.
Mixed dual mode models aren't implemented yet # todo
"""
function predict_concentration(dm::DualModeModel, pressures_mpa::AbstractVector)
    predictions = [predict_concentration(dm, pres_mpa) for pres_mpa in pressures_mpa]
    return predictions
end

"""
    predict_concentration(::AbstractVector{<:DualModeModel}, partial_pressures_mpa::AbstractVector{<:Number})
#todo 
Predict mixed gas concentrations using the dual mode mixing rule (langmuir-type competitive sorption) given a set of fit models and corresponding partial pressures.
e.g., 
    `predict_concentration([gas_1_model, gas_2_model, ...], [gas_1_pressure, gas_2_pressure, ...])`
!!! note
    `partial_pressures_mpa` is a bit of a misnomer. If the model was fit to fugacities, then fugacities should, of course, be specified (in MPa).
"""
function predict_concentration(dm::AbstractVector{<:DualModeModel}, partial_pressures_mpa::AbstractVector{<:Number})
    # todo
    throw(ErrorException("Not implemented yet, you should definitely complain to the devs about this."))
end


"""
Predict a pressure (**MPa**) given a concentration (**CC(STP)/CC(Polymer)**) according to the Dual Mode model.

The solved pressure for a single component Dual Mode model is:

``P = \\frac{\\sqrt{b^{2}\\left(c_h-x\\right)^{2}+2bk\\left(c+x\\right)+k^{2}}+b\\left(x-c\\right)-k}{2bk}``
"""
function predict_pressure(dm::DualModeModel, concentrations_cc_cc::Number)
    b = dm.b
    ch = dm.ch
    k = dm.kd
    c = concentrations_cc_cc
    
    if b * k != 0
        p = (sqrt((b^2)*(ch-c)^2+2*b*k*(ch + c) + k^2) + b*(c-ch) - k)/(2*b*k)
    elseif b == 0
        p = c / k
    else
        p = 0
    end
    
    return maximum(promote(0, p))
end

function MembraneBase.rss(dm::DualModeModel, isotherm::IsothermData)
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

function MembraneBase.rss(dm::DualModeModel, pres_or_fugs::AbstractVector, concentrations::AbstractVector)
    return rss(predict_concentration(dm, pres_or_fugs), concentrations)
end


function fit_dualmode_model(isotherm::IsothermData; uncertainty_method=nothing, use_fugacity=false, apply_weights=false)
    # see if isotherm is only a single component
    if isotherm.num_components != 1
        throw(ErrorException("The isotherm given has more than one component, this function only works for pure isotherms"))
    end

    if !apply_weights
        used_isotherm = strip_measurement_to_value(isotherm)
    else
        used_isotherm = isotherm
    end

    target = function(ch_b_kd)
        dmm = DualModeModel(ch_b_kd...; use_fugacity)
        err = rss(dmm, used_isotherm)
        if typeof(err) <: Measurement err = err.val end  # handle measurement types (we don't need them where we're going!)
        return err
    end
    lower = [0., 0., 0.]
    upper = [Inf, Inf, Inf]
    res = Optim.optimize(target, lower, upper, [0.5, 0.5, 0.5], Fminbox(BFGS()))
    
    if use_fugacity
        pressure_function = fugacities
        resampled_fitting_function = DualModeHelperFunctions.resampled_set_fitting_wrapper_fugacity
    else 
        pressure_function = partial_pressures
        resampled_fitting_function = DualModeHelperFunctions.resampled_set_fitting_wrapper_pressure
    end
    
    if !isnothing(uncertainty_method)
        ps = pressure_function(used_isotherm; component=1)
        cs = concentration(used_isotherm; component=1)
        data = collect(zip(ps, cs))
    end

    if uncertainty_method == :JackKnife
        corresponding_uncertainties = jackknife_uncertainty(resampled_fitting_function, data)
        uncertain_parameters = [Optim.minimizer(res)[i] ± corresponding_uncertainties[i] for i in 1:length(corresponding_uncertainties)]
        optimized_model = DualModeModel(uncertain_parameters...; use_fugacity)
    elseif uncertainty_method == :Bootstrap
        corresponding_uncertainties = bootstrap_uncertainty(resampled_fitting_function, data)  
        uncertain_parameters = [Optim.minimizer(res)[i] ± corresponding_uncertainties[i] for i in 1:length(corresponding_uncertainties)]
        optimized_model = DualModeModel(uncertain_parameters...; use_fugacity)  
    elseif isnothing(uncertainty_method)
        optimized_model = DualModeModel(Optim.minimizer(res)...; use_fugacity)
    else
        throw(ArgumentError("Invalid uncertainty_method: " * string(uncertainty_method)))
    end
    return optimized_model
end

"""
    fit_model(DualMode(), isotherm::IsothermData, [uncertainty_method=nothing], [apply_weights=false], [use_fugacity=false])
Fit the dual mode model to the pressures and concentrations present in the isotherm. 

Options
- For determining the uncertainty of the model parameters, the `:JackKnife`, and `:Bootstrap` methods are available. 
- `apply_weights` will use a weighted nonlinear regression method to solve the parameters, given that `Measurement` types are used somewhere in the data. 
- `use_fugacity` will fit the model to fugacities instead of pressures (they should be present in the isotherm data).  
"""
fit_model(::DualMode, isotherm::IsothermData; kwargs...) = fit_dualmode_model(isotherm; kwargs...)

function MembraneBase.strip_measurement_to_value(model::DualModeModel)
    return DualModeModel(
        strip_measurement_to_value(model.ch),
        strip_measurement_to_value(model.b),
        strip_measurement_to_value(model.kd),
        model.use_fugacity
    )
end

function thermo_factor(model::DualModeModel, pres_or_fug::Number, ρpol_g_cm3::Number, pen_mw::Number)
    ch = model.ch
    b = model.b 
    kd = model.kd
    p = pres_or_fug
    t1 = kd + ch * b / (1 + b * p)
    t2 = kd + ch * b / (1 + b * p)^2
    t3 = 1 + pen_mw / (ρpol_g_cm3 * MembraneBase.CC_PER_MOL_STP) * (t1 * p)
    return t1 / t2 * t3
end
