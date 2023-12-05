struct FloryHugginsDualMode end

struct FloryHugginsDualModeModel{CHIT, BT, KDT, PMVT} <: SorptionModel
    ch::CHIT  # units: CC/CC
    b::BT    # units: None
    chi::KDT  # units: None
    penetrant_molar_volume::PMVT # units: cm3/mol
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
    φ = flory_huggins_φ(fhdm, activity)
    flory_huggins_mode = concentration_from_estimated_volume_fraction(fhdm, φ) 
    langmuir_mode = fhdm.ch * fhdm.b * activity / (1 + fhdm.b * activity)
    return flory_huggins_mode + langmuir_mode
end

function estimate_volume_fraction(fhdm::FloryHugginsDualModeModel, concentration::Number)
    pen_volume = concentration * fhdm.penetrant_molar_volume / MembraneBase.CC_PER_MOL_STP  # ccstp/ccpol * cm3/mol / (ccstp/mol) = (cm3/mol * mol)/ccpol = cm3/ccpol = unitless
    vol_frac = pen_volume / (1 + pen_volume)  # unitless
    return vol_frac
end

function concentration_from_estimated_volume_fraction(fhdm::FloryHugginsDualModeModel, φ::number)
    φ / ((1-φ) * fhdm.penetrant_molar_volume / MembraneBase.CC_PER_MOL_STP)
end

flory_huggins_activity(fhdm::FloryHugginsDualModeModel, φ::Number) = exp(log(1-φ) + φ + fhdm.chi * φ^2)
flory_huggins_φ(fhdm::FloryHugginsDualModeModel, activity::Number) = Roots.find_zero((φ) -> activity - flory_huggins_activity(fhdm, φ), (0, 1))


function MembraneBase.rss(fhdm::FloryHugginsDualMode, isotherm::IsothermData)
    if isotherm.num_components == 1
        predictions = a_predict_concentration(fhdm, partial_pressures(isotherm, component=1))
        return MembraneBase.rss(concentration(isotherm, component=1), predictions)
    end
end

"""
    fit_model(FloryHugginsDualMode(), isotherm::IsothermData, [uncertainty_method=nothing], [apply_weights=false], [use_fugacity=false])
Fit the dual mode model to the pressures and concentrations present in the isotherm. 

Options
- For determining the uncertainty of the model parameters, the `:JackKnife`, and `:Bootstrap` methods are available. 
- `apply_weights` will use a weighted nonlinear regression method to solve the parameters, given that `Measurement` types are used somewhere in the data. 
- `use_fugacity` will fit the model to fugacities instead of pressures (they should be present in the isotherm data).  
"""
function fit_model(::FloryHugginsDualMode, isotherm::IsothermData, penetrant_partial_molar_volume; kwargs...) 
    # see if isotherm is only a single component
    if isotherm.num_components != 1
        throw(ErrorException("The isotherm given has more than one component, this function only works for pure isotherms"))
    end

    if !apply_weights
        used_isotherm = strip_measurement_to_value(isotherm)
    else
        used_isotherm = isotherm
    end

    target = function(ch_b_chi)
        fhdm = FloryHugginsDualModeModel(ch_b_chi..., penetrant_partial_molar_volume)
        err = rss(fhdm, used_isotherm)
        if typeof(err) <: Measurement err = err.val end  # handle measurement types (we don't need them where we're going!)
        return err
    end
    lower = [0., 0., 0.]
    upper = [Inf, Inf, Inf]
    res = Optim.optimize(target, lower, upper, [0.5, 0.5, 0.5], Fminbox(BFGS()))
    
    # if !isnothing(uncertainty_method)
    #     ps = pressure_function(used_isotherm; component=1)
    #     cs = concentration(used_isotherm; component=1)
    #     data = collect(zip(ps, cs))
    # end

    # if uncertainty_method == :JackKnife
    #     corresponding_uncertainties = jackknife_uncertainty(resampled_fitting_function, data)
    #     uncertain_parameters = [Optim.minimizer(res)[i] ± corresponding_uncertainties[i] for i in 1:length(corresponding_uncertainties)]
    #     optimized_model = DualModeModel(uncertain_parameters...; use_fugacity)
    # elseif uncertainty_method == :Bootstrap
    #     corresponding_uncertainties = bootstrap_uncertainty(resampled_fitting_function, data)  
    #     uncertain_parameters = [Optim.minimizer(res)[i] ± corresponding_uncertainties[i] for i in 1:length(corresponding_uncertainties)]
    #     optimized_model = DualModeModel(uncertain_parameters...; use_fugacity)  
    # elseif isnothing(uncertainty_method)
    optimized_model = DualModeModel(Optim.minimizer(res)...; use_fugacity)
    # else
        # throw(ArgumentError("Invalid uncertainty_method: " * string(uncertainty_method)))
    # end
    return optimized_model
end

function MembraneBase.strip_measurement_to_value(model::FloryHugginsDualMode)
    return FloryHugginsDualModeModel(
        strip_measurement_to_value(model.ch),
        strip_measurement_to_value(model.b),
        strip_measurement_to_value(model.chi),
        strip_measurement_to_value(model.penetrant_molar_volume)
    )
end