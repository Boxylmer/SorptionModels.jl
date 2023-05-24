
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
    
    sheet[wrp, wcp + 1] = "Pressure (MPa)"
    sheet[wrp, wcp + 2] = "Pressure σ"
    sheet[wrp, wcp + 3] = "Concentration (CC/CC)"
    sheet[wrp, wcp + 4] = "Concentration σ"
    sheet[wrp, wcp + 5] = "dp/dc (MPa/(CC/CC))"
    sheet[wrp, wcp + 6] = "dp/dc σ"
    sheet[wrp, wcp + 7] = "Frac. Dilation"
    sheet[wrp, wcp + 8] = "Frac. Dilation σ"
    sheet[wrp, wcp + 9] = "Continuous Dilation"
    sheet[wrp, wcp + 10] = "Cont. Dilation σ"   
    sheet[wrp, wcp + 11] = "d(ΔV/V0)/dp (1/MPa)"
    sheet[wrp, wcp + 12] = "d(ΔV/V0)/dp σ"
    sheet[wrp, wcp + 13] = "Partial molar volume (cm3/mol)"
    sheet[wrp, wcp + 14] = "Partial molar volume σ"

    wrp += 1
    
    sheet[wrp, wcp, dim=1] = collect(1:length(analysis.concentrations_cc_cc))
    
    write_vector_of_maybe_measurements(wrp, wcp + 1, analysis.pressures_mpa)
    write_vector_of_maybe_measurements(wrp, wcp + 3, analysis.concentrations_cc_cc)
    write_vector_of_maybe_measurements(wrp, wcp + 5, analysis.dp_dc)
    write_vector_of_maybe_measurements(wrp, wcp + 7, analysis.frac_dilations)
    write_vector_of_maybe_measurements(wrp, wcp + 9, analysis.continuous_dilations)
    write_vector_of_maybe_measurements(wrp, wcp + 11, analysis.dfracional_dilation_dp)
    write_vector_of_maybe_measurements(wrp, wcp + 13, analysis.partial_molar_volumes_cm3_mol)

end

function write_analysis(analysis::MolarVolumeAnalysis, filepath::AbstractString; name=DEFAULT_MOLAR_VOLUME_ANALYSIS_NAME)
    if !isfile(filepath)
        XLSX.openxlsx(filepath, mode="w") do _ end
    end
    
    XLSX.openxlsx(filepath, mode="rw") do xf 
        return write_analysis(analysis, xf; name=name)
    end
end
