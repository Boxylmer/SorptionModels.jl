struct PartialImmobilizationModel{HMDT, LMDT, RXVT, RYVT, PMT, PBT, TKT}
    henry_mode_diffusivity::HMDT
    langmuir_mode_diffusivity::LMDT
    regression_x_values::RXVT
    regression_y_values::RYVT

    pressures_mpa::PMT
    permeabilities_barrer::PBT
    temperature_k::TKT
    model::DualModeModel
end
"""
    PartialImmobilizationModel(::DualModeModel, ::AbstractVector{<:Number}, ::AbstractVector{<:Number}; [temperature_k])
Apply the method discussed in:

`P. Li, T.S. Chung, D.R. Paul, Gas sorption and permeation in PIM-1, Journal of Membrane Science. 432 (2013) 50–57. https://doi.org/10.1016/j.memsci.2013.01.009.`

to separate permeabilities into langmuir and henry mode diffusivities given a Dual Mode sorption model.



"""
function PartialImmobilizationModel(model::DualModeModel, pressures_mpa::AbstractVector{<:Number}, permeabilities_barrer::AbstractVector{<:Number}; temperature_k = nothing)
    if model.use_fugacity == false
        @warn "The Dual Mode model given was fit to pressures and not fugacities. The deconvolution will be less accurate as a result."
    end

    regression_x_values = model.ch .* model.b ./ (1 .+ model.b .* pressures_mpa) # cc/cc / MPa
    regression_y_values = permeabilities_barrer .* MembraneBase.CC_CC_MPA_CM2_S_PER_BARRER # barrer -> ((CC/CC)/MPa) * (cm^2/s) 
    
    slope, intercept = fit_linear_data(regression_x_values, regression_y_values) 
    # units: intercept: (CC/CC)/MPa * (cm^2/s); slope: cm^2/s
    
    henry_mode_diffusivity = intercept / model.kd
    langmuir_mode_diffusivity = slope

    return PartialImmobilizationModel(
        henry_mode_diffusivity, langmuir_mode_diffusivity,
        regression_x_values, regression_y_values,
        pressures_mpa, permeabilities_barrer, temperature_k, model
        )
end
