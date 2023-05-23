
DEFAULT_MOLAR_VOLUME_ANALYSIS_NAME = "Partial Molar Volume Analysis"

"""
write_analysis(analysis::MolarVolumeAnalysis, workbook::XLSX.XLSXFile; [name])
write_analysis(analysis::MolarVolumeAnalysis, filepath::AbstractString; [name])

Write a molar volume `analysis` out to a .xlsx `workbook`` or an .xlsx `filepath`. If a worksheet `name` is not provided, default names are created.

See [`MolarVolumeAnalysis`](@ref)
"""
function write_analysis(analysis::MolarVolumeAnalysis, workbook::XLSX.XLSXFile; name=DEFAULT_MOLAR_VOLUME_ANALYSIS_NAME)
    sheet = XLSX.addsheet!(workbook, name)
    wrp = 1  # writer row position
    wcp = 1  # writer col position

    function write_vector_of_maybe_measurements(row, col, vec; write_errors=true)
        sheet[row, col, dim=1] = strip_measurement_to_value(vec)
        if write_errors
            if eltype(vec) <: Measurement
                sheet[row, col+1, dim=1] = [vecval.err for vecval in vec]
            else
                sheet[row, col+1, dim=1] = ["No uncertainty" for _ in vec]
            end
        end
    end
    

    sheet[wrp, wcp] = "Step"
    sheet[wrp, wcp + 1] = "Concentration (CC/CC)"
    sheet[wrp, wcp + 2] = "Concentration σ"
    sheet[wrp, wcp + 3] = "dp/dc (MPa/(CC/CC))"
    sheet[wrp, wcp + 4] = "dp/dc σ"
    sheet[wrp, wcp + 5] = "d(ΔV/V0)/dp (1/MPa)"
    sheet[wrp, wcp + 6] = "d(ΔV/V0)/dp σ"
    sheet[wrp, wcp + 7] = "Partial molar volume (cm3/mol)"
    sheet[wrp, wcp + 8] = "Partial molar volume σ"

    wrp += 1
    
    sheet[wrp, wcp, dim=1] = collect(1:length(analysis.concentrations_cc_cc))
    write_vector_of_maybe_measurements(wrp, wcp + 1, analysis.concentrations_cc_cc)
    write_vector_of_maybe_measurements(wrp, wcp + 3, analysis.dp_dc)
    write_vector_of_maybe_measurements(wrp, wcp + 5, analysis.dfracional_dilation_dp)
    write_vector_of_maybe_measurements(wrp, wcp + 7, analysis.partial_molar_volumes_cm3_mol)

end

function write_analysis(analysis::MolarVolumeAnalysis, filepath::AbstractString; name=DEFAULT_MOLAR_VOLUME_ANALYSIS_NAME)
    if !isfile(filepath)
        XLSX.openxlsx(filepath, mode="w") do _ end
    end
    
    XLSX.openxlsx(filepath, mode="rw") do xf 
        return write_analysis(analysis, xf; name=name)
    end
end
