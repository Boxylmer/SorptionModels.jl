struct FloryHuggins end

struct FloryHugginsModel{KDT, PMVT} <: SorptionModel
    chi::KDT  # units: None
    penetrant_molar_volume::PMVT # units: cm3/mol
    activity_function::Union{Function, Missing}
end

function Base.show(io::IO, obj::FloryHugginsModel)
    kd = obj.chi
    capability = " with " * ((ismissing(obj.activity_function)) ? "no " : "" * "automatic activity conversion ability.")
    print(io, "Flory Huggins: χ = $kd" * capability)
end

estimate_volume_fraction(fh::FloryHugginsModel, concentration::Number) = flory_huggins_volume_fraction(fh.penetrant_molar_volume, concentration)

function flory_huggins_volume_fraction(penetrant_molar_volume::Number, concentration::Number)
    pen_volume = concentration * penetrant_molar_volume / MembraneBase.CC_PER_MOL_STP  # ccstp/ccpol * cm3/mol / (ccstp/mol) = (cm3/mol * mol)/ccpol = cm3/ccpol = unitless
    vol_frac = pen_volume / (1 + pen_volume)  # unitless
    return vol_frac
end

function flory_huggins_concentration_from_volume_fraction(fh::SorptionModel, φ::Number)
    φ / ((1-φ) * fh.penetrant_molar_volume / MembraneBase.CC_PER_MOL_STP)
end

flory_huggins_activity(fh::SorptionModel, φ::Number) = φ * exp((1-φ) + fh.chi * (1-φ)^2)
flory_huggins_φ(fh::SorptionModel, activity::Number) = Roots.find_zero((φ) -> activity - flory_huggins_activity(fh, φ), (0, 1))


function MembraneBase.rss(fh::FloryHugginsModel, activities::AbstractVector, concentrations::AbstractVector)
    predictions = [a_predict_concentration(fh, a) for a in activities]
    return MembraneBase.rss(concentrations, predictions)
    # # errs = ((concentration(isotherm, component=1) .-  predictions) ./ (concentration(isotherm, component=1) .+ 1e-16)).^2
    # errs = (concentration(isotherm, component=1) .- predictions).^2
    # return (sum(errs))
end



"""
"""
function fit_model(::FloryHuggins, activities::AbstractVector, penetrant_molar_volume::Number, concentrations::Base.AbstractVecOrTuple; activity_function::Function=missing, kwargs...) 
    target = function(chi)
        fh = FloryHugginsModel(chi..., penetrant_molar_volume)
        err = rss(fh, activities, concentrations) 
        if typeof(err) <: Measurement err = err.val end  # handle measurement types (we don't need them where we're going!)
        return err
    end
    # lower = [0., 0., -Inf]
    # upper = [Inf, Inf, Inf]
    # res = Optim.optimize(target, lower, upper, [1.0, 1.0, 0.], Fminbox(NelderMead()))

    res = Optim.optimize(target, [1.0], NelderMead())
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
    optimized_model = FloryHugginsModel(Optim.minimizer(res)..., penetrant_molar_volume, activity_function)
    # else
        # throw(ArgumentError("Invalid uncertainty_method: " * string(uncertainty_method)))
    # end
    return optimized_model
end

"""
"""
function fit_model(::FloryHuggins, pressures_mpa::Base.AbstractVecOrTuple, activity_function::Function, penetrant_molar_volume::Number, concentrations::Base.AbstractVecOrTuple; kwargs...) 
    return fit_model(
        FloryHuggins(), 
        activity_function.(pressures_mpa),
        penetrant_molar_volume,
        concentrations; 
        activity_function,
        kwargs...
    )
end

"""
    fit_model(FloryHuggins(), isotherm::IsothermData, [uncertainty_method=nothing], [apply_weights=false])
Fit the dual mode model to the pressures and concentrations present in the isotherm. 
If the isotherm contains activities, then they will be prioritized, otherwise it will look for a conversion function and use pressures.

Options
- For determining the uncertainty of the model parameters, the `:JackKnife`, and `:Bootstrap` methods are available. 
- `apply_weights` will use a weighted nonlinear regression method to solve the parameters, given that `Measurement` types are used somewhere in the data. 
- `use_fugacity` will fit the model to fugacities instead of pressures (they should be present in the isotherm data).  
"""
function fit_model(::FloryHuggins, isotherm::IsothermData, penetrant_molar_volume; apply_weights=false, activity_function::Function=missing, kwargs...) 
    # see if isotherm is only a single component
    if isotherm.num_components != 1
        throw(ErrorException("The isotherm given has more than one component, this function only works for pure isotherms"))
    end
    if !apply_weights
        used_isotherm = strip_measurement_to_value(isotherm)
    else
        used_isotherm = isotherm
    end

    if !isnothing(activities(isotherm))
        return fit_model(::FloryHuggins, activities(isotherm, component=1), penetrant_molar_volume, concentration(a, component=1); activity_function, kwargs...)
    else
        return fit_model(::FloryHuggins, partial_pressures(a, component=1), penetrant_molar_volume, concentration(a, component=1), activity_function; kwargs...)
    end
end


function MembraneBase.strip_measurement_to_value(model::FloryHuggins)
    return FloryHugginsModel(
        strip_measurement_to_value(model.chi),
        strip_measurement_to_value(model.penetrant_molar_volume),
        model.activity_function
    )
end