abstract type SorptionModel end

Base.broadcastable(x::SorptionModel) = Ref(x) # sorption models cannot be broadcasted

"""
    predict_concentration(model::SorptionModel, args...)

Predict concentration with a sorption model. Predictions are generally returned in the format of a vector of values corresponding to an input of a vector of pressures.
If the model typically uses activity, it must have an activity conversion function attached to it. 
"""
function predict_concentration(model::SorptionModel, potential::Number) end

"""
    a_predict_concentration(model::SorptionModel, args...)

Predict concentration with a sorption model that implements activity based methods. Predictions are generally returned in the format of a vector of values corresponding to an input of a vector of activities.
"""
function a_predict_concentration(model::SorptionModel, activity::Number) end



"""
    fit_model(model, ::IsothermData)

Fit a sorption model to an isotherm. If the data required for the model is not in the isotherm, an error message will be returned. 

Any keyword arguments get passed on to the specific model.
"""
function fit_model(_::SorptionModel, ::IsothermData) end

"Calculate the thermodynamic factor (see [ThermodynamicFactorAnalysis](@ref)) of a component in a polymer for a pure sorption model."
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
"""
    thermodynamic_factor(models::AbstractVector{<:SorptionModel}, pressures::AbstractVector{<:Number}, ρpol_g_cm3::Number, pen_mws::AbstractVector{<:Number}, zs::AbstractVector{<:Number}=ones(length(models))
(EXPERIMENTAL) Calculate the thermodynamic factor (see [ThermodynamicFactorAnalysis](@ref)) of a component in a polymer for a mixed sorption model.
This works by providing a vector of pure models, and implicitly uses their mixing rules to get back thermodynamic factors. Indices should represent components. 
"""
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

"""
    dc_dp(model::SorptionModel, p_mpa::Number)
Calculate the derivative of concentration (cc/cc) with respect to pressure (MPa).
"""
function dc_dp(model::SorptionModel, p_mpa::Number)
    return ForwardDiff.derivative.(p -> predict_concentration(model, p), p_mpa) # mpa / cc/cc
end

"""
    dp_dc(model::SorptionModel, p_mpa::Number)
Calculate the derivative of pressure (MPa) with respect to concentration (cc/cc).
"""
function dp_dc(model::SorptionModel, p_mpa::Number)
    return 1 / dc_dp(model, p_mpa) # mpa / cc/cc
end 