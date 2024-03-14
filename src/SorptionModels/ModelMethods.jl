abstract type SorptionModel end

"""
    predict_concentration(model::SorptionModel, args...)

Predict concentration with a sorption model. Predictions are generally returned in the format of a vector of values corresponding to an input of a vector of pressures.
"""
function predict_concentration(model::SorptionModel, potential::Number) end

"""
    fit_model(model, ::IsothermData)

Fit a sorption model to an isotherm. If the data required for the model is not in the isotherm, an error message will be returned. 

Any keyword arguments get passed on to the specific model.
"""
function fit_model(_::SorptionModel, ::IsothermData) end

function thermodynamic_factor(model::SorptionModel, pressure::Number, ρpol_g_cm3::Number, pen_mw::Number, z::Number=1.0)

    ln_w(pressure_mpa) = log(
        ccpen_per_ccpol_to_mass_fractions(
            predict_concentration(model, pressure_mpa), ρpol_g_cm3, [pen_mw]
            )[2]
        ) # verified

    ln_c(pressure_mpa) = log(predict_concentration(model, pressure_mpa)) # verified

    c = predict_concentration(model, pressure) # verified
    dlnw_dp = ForwardDiff.derivative(ln_w, pressure) # verified
    dlnc_dp = ForwardDiff.derivative(ln_c, pressure) # verified
    dc_dp = ForwardDiff.derivative(p -> predict_concentration(model, p), pressure) # verified

    # α = z * d(ln(ω))/d(ln(c)) * (c/p) * dp/dc
    # α = z * d(ln(ω))/dp / d(ln(c))/dp * (c/p) * 1 / (dc/dp)
    # @show dlnw_dp ./ dlnc_dp # verified

    α = z * dlnw_dp / dlnc_dp * (c / pressure) * (1 / dc_dp) 
    return α
end