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
        )
    dlnw_dp = ForwardDiff.derivative(ln_w, pressure) 

    # α = z * dlnc_dp / dlnw_dp * (c / pressure) * (1 / dc_dp) = (Z/p)*(dp/dln(w))
    α = z / pressure * 1 / dlnw_dp
    return α
end

# mixed gas case, requires the SorptionModel supposed mixing 
function thermodynamic_factor(
    models::AbstractVector{<:SorptionModel}, 
    pressures::AbstractVector{<:Number}, 
    ρpol_g_cm3::Number, 
    pen_mws::AbstractVector{<:Number}, 
    zs::AbstractVector{<:Number}=ones(length(models))
    )

    dlnw_dps = Vector(undef, length(models)) 
    
    working_pressure_vector = Vector{Number}(undef, length(models)) 
    working_pressure_vector .= pressures
    for i in eachindex(models)
        function ln_w(p)
            temp = working_pressure_vector[i]
            working_pressure_vector[i] = p
            lnw = log(
            ccpen_per_ccpol_to_mass_fractions(
                predict_concentration(models, working_pressure_vector), ρpol_g_cm3, pen_mws
                )[1 + i]
            )
            working_pressure_vector[i] = temp
            return lnw
        end
        dlnw_dp = ForwardDiff.derivative(ln_w, pressures[i]) 
        dlnw_dps[i] = dlnw_dp
    end

    # α = z * dlnc_dp / dlnw_dp * (c / pressure) * (1 / dc_dp) = (Z/p)*(dp/dln(w))
    αs = zs ./ pressures .* 1 ./ dlnw_dps
    return αs
end