struct FloryHugginsDualMode end

nparams(::FloryHugginsDualMode) = 3

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

# """
#     infinite_dilution_solubility(dm::FloryHugginsDualModeModel)
# Get the infinite dilution solubility coefficient in (**(CC/CC)**) (activity is unitless)
# """
# infinite_dilution_solubility(dm::FloryHugginsDualModeModel) = nothing # TODO


"""
    a_predict_concentration(::FloryHugginsDualModeModel, activity::Number)
Predict the concentration of a penetrant given a flory huggins dual mode model and activity.
"""
function a_predict_concentration(fhdm::FloryHugginsDualModeModel, activity::Number)
    φ = flory_huggins_φ(fhdm, activity)
    flory_huggins_mode = flory_huggins_concentration_from_volume_fraction(fhdm, φ) 
    langmuir_mode = fhdm.ch * fhdm.b * activity / (1 + fhdm.b * activity)
    return flory_huggins_mode + langmuir_mode
end

estimate_volume_fraction(fhdm::FloryHugginsDualModeModel, concentration::Number) = flory_huggins_volume_fraction(fhdm.penetrant_molar_volume, concentration)

# TODO this needs to be generalized
function MembraneBase.rss(fhdm::FloryHugginsDualModeModel, isotherm::IsothermData)
    predictions = [a_predict_concentration(fhdm, a) for a in activities(isotherm, component=1)]
    return MembraneBase.rss(concentration(isotherm, component=1), predictions)
    # # errs = ((concentration(isotherm, component=1) .-  predictions) ./ (concentration(isotherm, component=1) .+ 1e-16)).^2
    # errs = (concentration(isotherm, component=1) .- predictions).^2
    # return (sum(errs))
end

"""
    fit_model(FloryHugginsDualMode(), isotherm::IsothermData, [uncertainty_method=nothing], [apply_weights=false])
Fit the dual mode model to the pressures and concentrations present in the isotherm. 

Options
- For determining the uncertainty of the model parameters, the `:JackKnife`, and `:Bootstrap` methods are available. 
- `apply_weights` will use a weighted nonlinear regression method to solve the parameters, given that `Measurement` types are used somewhere in the data. 
- `use_fugacity` will fit the model to fugacities instead of pressures (they should be present in the isotherm data).  
"""
function fit_model(::FloryHugginsDualMode, isotherm::IsothermData, penetrant_molar_volume; apply_weights=false, kwargs...) 
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
        fhdm = FloryHugginsDualModeModel(ch_b_chi..., penetrant_molar_volume)
        err = rss(fhdm, used_isotherm)
        if typeof(err) <: Measurement err = err.val end  # handle measurement types (we don't need them where we're going!)
        return err
    end
    # lower = [0., 0., -Inf]
    # upper = [Inf, Inf, Inf]
    # res = Optim.optimize(target, lower, upper, [1.0, 1.0, 0.], Fminbox(NelderMead()))

    res = Optim.optimize(target, [1.0, 1.0, 1.0], NelderMead())
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
    optimized_model = FloryHugginsDualModeModel(Optim.minimizer(res)..., penetrant_molar_volume)
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