struct DualModeDesorption{SMT, DMT}
    sorbing_model::SMT
    desorbing_model::DMT
end

function DualModeDesorption(isotherm::IsothermData; use_fugacity=false)
    sorbing_isotherm = increasing_concentration(isotherm)
    desorbing_isotherm = remove_increasing_concentration_steps(isotherm)

    

    if num_steps(sorbing_isotherm) <= 1 || num_steps(desorbing_isotherm) <= 1
        throw(BoundsError("Isotherm either did not contain any sorption steps or did not contain any desorption steps."))
    end

    function target(params)
        sorption_model = DualModeModel(params[1], params[2], params[3]; use_fugacity)
        desorption_model = DualModeModel(params[4], params[2], params[3]; use_fugacity)
        sorption_rss = rss(sorption_model, sorbing_isotherm)
        desorption_rss = rss(desorption_model, desorbing_isotherm)
        return sorption_rss + desorption_rss
    end

    lower = [0., 0., 0.]
    upper = [Inf, Inf, Inf]
    res = Optim.optimize(target, lower, upper, [0.5, 0.5, 0.5], Fminbox(BFGS()))
    optimized_params = Optim.minimizer(res)

end
