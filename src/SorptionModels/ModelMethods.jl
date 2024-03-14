abstract type SorptionModel end

"""
    predict_concentration(model::SorptionModel, args...)

Predict concentration with a sorption model. Predictions are generally returned in the format of a vector of values corresponding to an input of a vector of pressures.
"""
function predict_concentration end

"""
    fit_model(model, ::IsothermData)

Fit a sorption model to an isotherm. If the data required for the model is not in the isotherm, an error message will be returned. 

Any keyword arguments get passed on to the specific model.
"""
function fit_model(model::SorptionModel, ::IsothermData) end



function thermodynamic_factor(model::SorptionModels, pressure::Number, ρpol_g_cm3::Number, pen_mw::Number, z::Number=1.0)
    # # looking for ∂ln(a)/∂ln(ω)
    ln_a(pressure_mpa) = log(activity_function(pressure_mpa))
    ln_w(pressure_mpa) = log(
        ccpen_per_ccpol_to_mass_fractions(
            predict_concentration(model, pressure_mpa), ρpol_g_cm3, [pen_mw]
        )[2]
    )

    ln_c(pressure_mpa) = log(predict_concentration(model, pressure_mpa))
    c = predict_concentration(model, pressure_mpa)

    dlnc_dp = ForwardDiff.derivative(ln_w, pressure)
    dlna_dp = ForwardDiff.derivative(ln_a, pressure)
    dlnw_dp = ForwardDiff.derivative(ln_w, pressure)
    
    # α = z * d(ln(ω))/d(ln(c)) * (c/p) * dp/dc
    # α = z * d(ln(ω))/dp / d(ln(c))/dp * (c/p) * 1 / (dc/dp)
    α = dlna_dp ./ dlnw_dp * z

    return α
end
end