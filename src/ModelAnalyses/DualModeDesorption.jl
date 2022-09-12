struct DualModeDesorption{SMT, DMT}
    sorbing_model::SMT
    desorbing_model::DMT
    max_p::Float64
    isotherm::IsothermData
end

function DualModeDesorption(isotherm::IsothermData; use_fugacity=false, uncertainty_method=nothing, naive=false, share_b=true, verbose=false)
    sorbing_isotherm = increasing_concentration(isotherm)
    desorbing_isotherm = remove_increasing_concentration_steps(isotherm)

    if use_fugacity
        max_p = maximum(fugacity(isotherm; component=1))
    else
        max_p = maximum(partial_pressures(isotherm; component=1))
    end
    
    if num_steps(sorbing_isotherm) <= 1 || num_steps(desorbing_isotherm) <= 1
        if verbose
            @warn "Isotherm either did not contain any sorption steps or did not contain any desorption steps."
        end
        return nothing
    end

    if naive
        return DualModeDesorption(
            fit_model(DualMode(), sorbing_isotherm; use_fugacity, uncertainty_method), 
            fit_model(DualMode(), desorbing_isotherm; use_fugacity, uncertainty_method), 
            measurement(max_p).val,
            isotherm
        )
    end

    concs = concentration(isotherm; component=1)

    if use_fugacity
        pres = fugacity(isotherm; component=1)
    else
        pres = partial_pressures(isotherm; component=1)
    end

    data = collect(zip(pres, concs))

    fitter = make_dualmode_desorption_fitting_function(use_fugacity, share_b)
    params = fitter(data)

    if uncertainty_method == :JackKnife
        uncertainties = jackknife_uncertainty(fitter, data)
        params = params .± uncertainties
    elseif isnothing(uncertainty_method)
        
    else
        throw(ArgumentError("Uncertainty method $uncertainty_method is invalid. This method supports :JackKnife or nothing."))
    end

    sorbing_model = DualModeModel(params[1], params[2], params[3]; use_fugacity)
    desorbing_model = DualModeModel(params[4], params[2], params[5]; use_fugacity)
    return DualModeDesorption(sorbing_model, desorbing_model, measurement(max_p).val, isotherm)

end

function make_dualmode_desorption_fitting_function(use_fugacity, share_b=true)::Function
    function fitting_function(param_tuples)::Vector{Float64}
        # convert the data tuples
        pressures, concs = zip(param_tuples...)

        # construct the isotherm
        if use_fugacity
            iso = IsothermData(fugacities_mpa = pressures, concentrations_cc = concs)
        else
            iso = IsothermData(partial_pressures_mpa = pressures, concentrations_cc = concs)
        end

        # split the isotherm into a sorption and desorption component
        sorbing_isotherm = increasing_concentration(iso)
        desorbing_isotherm = remove_increasing_concentration_steps(iso)
        if num_steps(sorbing_isotherm) <= 1 || num_steps(desorbing_isotherm) <= 1
            throw(BoundsError("Isotherm either did not contain any sorption steps or did not contain any desorption steps."))
        end
        
    
        # optimize the target
        if share_b
            function target(params)
                sorption_model = DualModeModel(params[1], params[2], params[3]; use_fugacity)
                desorption_model = DualModeModel(params[4], params[2], params[5]; use_fugacity)
                sorption_rss = rss(sorption_model, sorbing_isotherm)
                desorption_rss = rss(desorption_model, desorbing_isotherm)
                return sorption_rss + desorption_rss
            end

            lower = [0., 0., 0., 0., 0.]
            upper = [Inf, Inf, Inf, Inf, Inf]
            init = [0.5, 0.5, 0.5, 0.5, 0.5]
            res = Optim.optimize(target, lower, upper, init, Fminbox(BFGS()))
            optimized_params = Optim.minimizer(res)
            return optimized_params
        else
            m = fit_model(DualMode(), sorbing_isotherm; use_fugacity)

            function static_b_target(params)
                desorption_model = DualModeModel(params[1], m.b, params[2]; use_fugacity)
                desorption_rss = rss(desorption_model, desorbing_isotherm)
                return desorption_rss
            end

            lower = [0., 0.]
            upper = [Inf, Inf]
            init = [0.5, 0.5]
            res = Optim.optimize(static_b_target, lower, upper, init, Fminbox(BFGS()))
            optimized_params = Optim.minimizer(res)
            return [m.ch, m.b, m.kd, optimized_params[1], optimized_params[2]]
        end
    end

    return fitting_function
end

function predict_concentration(model::DualModeDesorption, isotherm::IsothermData)
    sorbing_isotherm = increasing_concentration(isotherm)
    desorbing_isotherm = remove_increasing_concentration_steps(isotherm)
    if model.sorbing_model.use_fugacity && model.desorbing_model.use_fugacity
        sps = fugacities(sorbing_isotherm; component=1)
        dps = fugacities(desorbing_isotherm; component=1)
    elseif !model.sorbing_model.use_fugacity && !model.desorbing_model.use_fugacity 
        sps = partial_pressures(sorbing_isotherm; component=1)
        dps = partial_pressures(desorbing_isotherm; component=1)
    else
        throw(ErrorException("The Dual Mode Desorption Analysis was not valid. One model used fugacity and another did not. Was this analysis manually created?"))
    end
    return predict_concentration(model, sps, dps)
end

function predict_concentration(model::DualModeDesorption, sorbing_pressures, desorbing_pressures)
    return predict_concentration(model.sorbing_model, sorbing_pressures), predict_concentration(model.desorbing_model, desorbing_pressures)
end

function predict_concentration(model::DualModeDesorption, pressures)
    return predict_concentration(model.sorbing_model, pressures), predict_concentration(model.desorbing_model, pressures)
end

function Base.show(io::IO, obj::DualModeDesorption)
    print(io, "Sorbing: " * string(obj.sorbing_model) * ", Desorbing: " * string(obj.desorbing_model))
end

using_fugacity(dmd::DualModeDesorption) = dmd.sorbing_model.use_fugacity