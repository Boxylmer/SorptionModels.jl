DEFAULT_PARTIAL_IMMOBILIZATION_ANALYSIS_NAME = "Partial Immobilization Analysis"

"""
    write_analysis(analysis::PartialImmobilizationModel, workbook::XLSX.XLSXFile; [name])
    write_analysis(analysis::PartialImmobilizationModel, filepath::AbstractString; [name])

Write a partial immobilization `analysis` out to a .xlsx `workbook`` or an .xlsx `filepath`. If a worksheet `name` is not provided, default names are created.

See [`PartialImmobilizationModel`](@ref)
"""
function write_analysis(analysis::PartialImmobilizationModel, workbook::XLSX.XLSXFile; name=DEFAULT_PARTIAL_IMMOBILIZATION_ANALYSIS_NAME)
    sheet = XLSX.addsheet!(workbook, name)
    wrp = 1  # writer row position
    wcp = 1  # writer col position

    function write_vector_of_maybe_measurements(row, col, vec; write_errors=true)
        sheet[row, col, dim=1] = strip_measurement_to_value(vec)
        if write_errors
            if eltype(vec) <: Measurement
                sheet[row, col+1, dim=1] = [vecval.err for vecval in vec]
            else
                sheet[row, col+1, dim=1] = ["No uncertainty"]
            end
        end
    end

    function write_maybe_measurement(row, col, meas)
        sheet[row, col] = strip_measurement_to_value(meas)
        if typeof(meas) <: Measurement
            sheet[row, col+1] = meas.err 
        else
            sheet[row, col+1] = "No uncertainty"
        end
    end
    
    sheet[wrp, wcp] = "Temperature (K)"
    if analysis.model.use_fugacity == false
        sheet[wrp, wcp + 1] = "P (MPa)"
        sheet[wrp, wcp + 2] = "P σ"
    else
        sheet[wrp, wcp + 1] = "Fugacity (MPa)"
        sheet[wrp, wcp + 2] = "Fugacity σ"
    end
    
    sheet[wrp, wcp + 3] = "C (CC/CC)"
    sheet[wrp, wcp + 4] = "C σ"

    sheet[wrp, wcp + 5] = "Perm. (Barrer)"
    sheet[wrp, wcp + 6] = "Perm. σ"
    sheet[wrp, wcp + 7] = "DH (cm^2/s)"
    sheet[wrp, wcp + 8] = "DH σ"   
    sheet[wrp, wcp + 9] = "DD (cm^2/s)"
    sheet[wrp, wcp + 10] = "DD σ"   
    sheet[wrp, wcp + 11] = "F (Dh/Dd)"
    sheet[wrp, wcp + 12] = "F σ"

    if analysis.model.use_fugacity == false
        sheet[wrp, wcp + 13] = "ch'*b/(1+bp)"
        sheet[wrp, wcp + 14] = "ch'*b/(1+bp) σ"
    else
        sheet[wrp, wcp + 13] = "ch'*b/(1+bf)"
        sheet[wrp, wcp + 14] = "ch'*b/(1+bf) σ"
    end

    wrp += 1

    if !isnothing(analysis.temperature_k)
        sheet[wrp, wcp] = typeof(analysis.temperature_k) <: Measurement ? analysis.temperature_k.val : analysis.temperature_k
    else
        sheet[wrp, wcp] = "Not specified"
    end
    write_vector_of_maybe_measurements(wrp, wcp + 1, analysis.pressures_mpa)
    write_vector_of_maybe_measurements(wrp, wcp + 3, predict_concentration(analysis.model, analysis.pressures_mpa))
    write_maybe_measurement(wrp, wcp + 5, analysis.permeabilities_barrer)
    write_maybe_measurement(wrp, wcp + 7, analysis.langmuir_mode_diffusivity)
    write_maybe_measurement(wrp, wcp + 9, analysis.henry_mode_diffusivity)
    write_maybe_measurement(wrp, wcp + 11, analysis.f)
    write_vector_of_maybe_measurements(wrp, wcp + 13, analysis.regression_x_values)
end

function write_analysis(analysis::PartialImmobilizationModel, filepath::AbstractString; name=DEFAULT_PARTIAL_IMMOBILIZATION_ANALYSIS_NAME)
    if !isfile(filepath)
        XLSX.openxlsx(filepath, mode="w") do _ end
    end
    
    XLSX.openxlsx(filepath, mode="rw") do xf 
        return write_analysis(analysis, xf; name=name)
    end
end
