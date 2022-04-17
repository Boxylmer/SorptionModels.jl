"""
write_analysis(analysis::PartialImmobilizationModel, filepath::AbstractString; [analysis_name])
write_dual_mode_diffusivity_deconvolution(analyses::AbstractVector{PartialImmobilizationModel}, filepath::AbstractString; [analysis_names])

Write a Partial Immobilization Model `analysis` out to a .csv titled by `filepath`. If analysis names are not provided, default names are created.
Methods are provided for a single analysis, or vectors of analyses (all written to one file).

See [`PartialImmobilizationModel`](@ref)
"""
function write_analysis(
    analysis::PartialImmobilizationModel, 
    filepath::AbstractString; 
    analyses_name = missing)
    write_analysis([analysis], filepath; analyses_names = [analyses_name])
end

function write_analysis(
    analyses::AbstractVector{<:PartialImmobilizationModel},
    filepath::AbstractString;
    analyses_names::AbstractVector = missing)
    permeability_deconvolution_results = []

    if ismissing(analyses_names)
        analyses_names = "Analysis " .* string.(1:length(analyses))
    end

    for (idx, analysis) in enumerate(analyses)
        langmuir_diffusivity, henry_diffusivity = analysis.langmuir_mode_diffusivity, analysis.henry_mode_diffusivity
        f = langmuir_diffusivity / henry_diffusivity
        x_values = [measurement(xval).val for xval in analysis.regression_x_values]

        # ch_b_over_1_plus_fp = model.ch * model.b / (1 .+ model.b .* corresponding_fugacities)
        push!(permeability_deconvolution_results, [analyses_names[idx]])
        push!(permeability_deconvolution_results, ["Temperature (K): ", strip_measurement_to_value(analysis.temperature_k)])
        if analysis.model.use_fugacity == false
            push!(permeability_deconvolution_results, ["Pressures (MPa): ", analysis.pressures_mpa...])
        else
            push!(permeability_deconvolution_results, ["Fugacities (MPa): ", analysis.pressures_mpa...])
        end
        push!(permeability_deconvolution_results, ["Permeabilities (Barrer): ", analysis.permeabilities_barrer...])
        push!(permeability_deconvolution_results, ["Langmuir Diffusivity (cm^2/s): ", measurement(langmuir_diffusivity).val, "+/-", measurement(langmuir_diffusivity).err])
        push!(permeability_deconvolution_results, ["Henry Diffusivity (cm^2/s): ", measurement(henry_diffusivity).val, "+/-", measurement(henry_diffusivity).err])
        push!(permeability_deconvolution_results, ["F (Dh/Dd): ", measurement(f).val, "+/-", measurement(f).err])
        push!(permeability_deconvolution_results, ["ch'*b/(1+bf)", x_values...])
    end
    writedlm(
        joinpath(
            @__DIR__, 
            filepath, 
        ), permeability_deconvolution_results, ",")

end