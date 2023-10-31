"""
    write_analysis(analysis::VantHoffDualModeAnalysis, workbook::XLSX.XLSXFile; [name])
    write_analysis(analysis::VantHoffDualModeAnalysis, filepath::AbstractString; [name])

Write a Vant Hoff Dual Mode `analysis` out to a .xlsx `workbook`` or an .xlsx `filepath`. If a worksheet `name` is not provided, default names are created.

See [`VantHoffDualModeAnalysis`](@ref)
"""
function write_analysis(analysis::VantHoffDualModeAnalysis, workbook::XLSX.XLSXFile; name="Van't Hoff Dual Mode Analysis")
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
    function write_value_of_maybe_measurements(row, col, val; write_errors=true)
        sheet[row, col] = strip_measurement_to_value(val)
        if write_errors
            if typeof(val) <: Measurement
                sheet[row, col+1] = val.err
            else
                sheet[row, col+1] = "No uncertainty"
            end
        end
    end
    function write_parameter_headers(row_start, col_start; row_col_iter=(1, 0)) 
        rowpos = 0
        colpos = 0
        rowpos += row_col_iter[1]; colpos += row_col_iter[2]
        sheet[row_start + rowpos, col_start + colpos] = "CH'_0 (CC/CC)"
        rowpos += row_col_iter[1]; colpos += row_col_iter[2]
        sheet[row_start + rowpos, col_start + colpos] = "b_0 (1/MPa)"
        rowpos += row_col_iter[1]; colpos += row_col_iter[2]
        sheet[row_start + rowpos, col_start + colpos] = "kd_0 ((CC/CC)/MPa)"
        rowpos += row_col_iter[1]; colpos += row_col_iter[2]
        sheet[row_start + rowpos, col_start + colpos] = "m_CH` (1/K)"
        rowpos += row_col_iter[1]; colpos += row_col_iter[2]
        sheet[row_start + rowpos, col_start + colpos] = "ΔHb (J/mol)"
        rowpos += row_col_iter[1]; colpos += row_col_iter[2]
        sheet[row_start + rowpos, col_start + colpos] = "ΔHkd (J/mol)"
    end
    function write_matrix_at_position(row, col, mat)
        for idx in CartesianIndices(mat)
            sheet[row + idx[1] - 1, col + idx[2] - 1] = mat[idx[1], idx[2]]
        end
    end

    sheet[wrp, wcp] = "Parameters"
    sheet[wrp, wcp + 1] = "Initial"
    sheet[wrp, wcp + 2] = "Final"
    sheet[wrp, wcp + 3] = "Final err"

    write_parameter_headers(wrp, wcp)

    write_value_of_maybe_measurements(wrp+1, wcp+1, analysis.ch0_initial; write_errors=false); write_value_of_maybe_measurements(wrp+1, wcp+2, analysis.ch0_final)
    write_value_of_maybe_measurements(wrp+2, wcp+1, analysis.b0_initial; write_errors=false); write_value_of_maybe_measurements(wrp+2, wcp+2, analysis.b0_final)
    write_value_of_maybe_measurements(wrp+3, wcp+1, analysis.kd0_initial; write_errors=false); write_value_of_maybe_measurements(wrp+3, wcp+2, analysis.kd0_final)
    write_value_of_maybe_measurements(wrp+4, wcp+1, analysis.mch_initial; write_errors=false); write_value_of_maybe_measurements(wrp+4, wcp+2, analysis.mch_final)
    write_value_of_maybe_measurements(wrp+5, wcp+1, analysis.ΔHb_initial; write_errors=false); write_value_of_maybe_measurements(wrp+5, wcp+2, analysis.ΔHb_final)
    write_value_of_maybe_measurements(wrp+6, wcp+1, analysis.ΔHkd_initial; write_errors=false); write_value_of_maybe_measurements(wrp+6, wcp+2, analysis.ΔHkd_final)

    sheet[wrp+7, wcp] = "RSS"
    write_value_of_maybe_measurements(wrp+7, wcp+1, analysis.initial_rss; write_errors=false); write_value_of_maybe_measurements(wrp+7, wcp+2, analysis.final_rss)

    wcp+=5
    sheet[wrp, wcp] = "Covariance Matrix"
    write_parameter_headers(wrp, wcp)
    write_parameter_headers(wrp, wcp; row_col_iter=(0, 1))
    write_matrix_at_position(wrp + 1, wcp + 1, analysis.covariance_matrix)

    wcp+=9
    sheet[wrp, wcp] = "1/T"
    write_vector_of_maybe_measurements(wrp + 1, wcp, 1 ./ [temperature(iso) for iso in analysis.given_isotherms]; write_errors=false)
    sheet[wrp, wcp + 1] = "ch (standalone)"
    write_vector_of_maybe_measurements(wrp + 1, wcp + 1, [model.ch for model in analysis.standalone_models]; write_errors=false)
    sheet[wrp, wcp + 2] = "ln(b) (standalone)"
    write_vector_of_maybe_measurements(wrp + 1, wcp + 2, log.([model.b for model in analysis.standalone_models]); write_errors=false)
    sheet[wrp, wcp + 3] = "ln(kd) (standalone)"
    write_vector_of_maybe_measurements(wrp + 1, wcp + 3, log.([model.kd for model in analysis.standalone_models]); write_errors=false)
    sheet[wrp, wcp + 4] = "ch (initial)"
    write_vector_of_maybe_measurements(wrp + 1, wcp + 4, [model.ch for model in analysis.initial_models]; write_errors=false)
    sheet[wrp, wcp + 5] = "ln(b) (initial)"
    write_vector_of_maybe_measurements(wrp + 1, wcp + 5, log.([model.b for model in analysis.initial_models]); write_errors=false)
    sheet[wrp, wcp + 6] = "ln(kd) (initial)"
    write_vector_of_maybe_measurements(wrp + 1, wcp + 6, log.([model.kd for model in analysis.initial_models]); write_errors=false)
    sheet[wrp, wcp + 7] = "ch (final)"
    write_vector_of_maybe_measurements(wrp + 1, wcp + 7, [model.ch for model in analysis.final_models]; write_errors=false)
    sheet[wrp, wcp + 8] = "ln(b) (final)"
    write_vector_of_maybe_measurements(wrp + 1, wcp + 8, log.([model.b for model in analysis.final_models]); write_errors=false)
    sheet[wrp, wcp + 9] = "ln(kd) (final)"
    write_vector_of_maybe_measurements(wrp + 1, wcp + 9, log.([model.kd for model in analysis.final_models]); write_errors=false)
    wcp=1

    wrp+=7

    wrp+= max(2, length(analysis.given_isotherms) - 8)  # move the writer by the longest vector
    sheet[wrp, wcp] = "Isotherms"
    wrp+=1


    for idx in eachindex(analysis.final_models, analysis.initial_models, analysis.final_predicted_isotherms, analysis.initial_predicted_isotherms)
        sheet[wrp, wcp] = "Temp (K):"
        write_value_of_maybe_measurements(wrp, wcp+1, temperature(analysis.given_isotherms[idx]))
        
        # isotherm stuff
        wrp+=1
        if analysis.use_fugacity sheet[wrp, wcp] = "Fugacity (MPa)"
        else sheet[wrp, wcp] = "Pressure (MPa)" end
        sheet[wrp, wcp + 1], sheet[wrp, wcp + 2], sheet[wrp, wcp + 3], sheet[wrp, wcp + 4] = "Err (MPa)", "CC/CC", "CC/CC err", "DM standalone (CC/CC)"
        sheet[wrp, wcp + 5], sheet[wrp, wcp + 6], sheet[wrp, wcp + 7], sheet[wrp, wcp + 8], sheet[wrp, wcp + 9]  = "DM init (CC/CC)", "DM final (CC/CC)", "err (CC/CC)", "DM Henry Contrib (CC/CC)", "err (CC/CC)"
        sheet[wrp, wcp + 10], sheet[wrp, wcp + 11] = "DM Langmuir Contrib (CC/CC)", "err (CC/CC)"
        
        
        wrp+=1
        
        if analysis.use_fugacity
            used_pressures = fugacities(analysis.given_isotherms[idx]; component=1)
        else
            used_pressures = partial_pressures(analysis.given_isotherms[idx]; component=1)
        end
        write_vector_of_maybe_measurements(wrp, wcp, fugacities(analysis.given_isotherms[idx]; component=1))
        write_vector_of_maybe_measurements(wrp, wcp + 2, concentration(analysis.given_isotherms[idx]; component=1))
        write_vector_of_maybe_measurements(wrp, wcp + 4, strip_measurement_to_value(predict_concentration(analysis.standalone_models[idx], used_pressures)))
        write_vector_of_maybe_measurements(wrp, wcp + 5, strip_measurement_to_value(predict_concentration(analysis.initial_models[idx], used_pressures)))
        write_vector_of_maybe_measurements(wrp, wcp + 6, predict_concentration(analysis.final_models[idx], used_pressures))
        write_vector_of_maybe_measurements(wrp, wcp + 8, henry_mode_concentration(analysis.final_models[idx], used_pressures))
        write_vector_of_maybe_measurements(wrp, wcp + 10, langmuir_mode_concentration(analysis.final_models[idx], used_pressures))
        
        # dual mode fittings
        wrp+=0
        wcp+=14
        sheet[wrp - 1, wcp] = "Dual Mode Parameters"
        sheet[wrp + 1, wcp], sheet[wrp + 2, wcp], sheet[wrp + 3, wcp] = "CH' (CC/CC)", "b (1/MPa)", "kd ((CC/CC)/MPa)"
        sheet[wrp, wcp + 1], sheet[wrp, wcp + 2], sheet[wrp, wcp + 3], sheet[wrp, wcp + 4], sheet[wrp, wcp + 5] = "standalone", "standalone err", "initial", "final", "final err"
        wrp+=1
        write_vector_of_maybe_measurements(wrp, wcp + 1, [analysis.standalone_models[idx].ch, analysis.standalone_models[idx].b, analysis.standalone_models[idx].kd])
        write_vector_of_maybe_measurements(wrp, wcp + 3, [analysis.initial_models[idx].ch, analysis.initial_models[idx].b, analysis.initial_models[idx].kd]; write_errors=false)
        write_vector_of_maybe_measurements(wrp, wcp + 4, [analysis.final_models[idx].ch, analysis.final_models[idx].b, analysis.final_models[idx].kd])

        wrp+= 1 + num_steps(analysis.given_isotherms[idx])
        wcp = 1

    end

end

function write_analysis(analysis::VantHoffDualModeAnalysis, filepath::AbstractString; name="Van't Hoff Dual Mode Analysis")
    if !isfile(filepath)
        XLSX.openxlsx(filepath, mode="w") do _ end
    end
    
    XLSX.openxlsx(filepath, mode="rw") do xf 
        return write_analysis(analysis, xf; name=name)
    end
end
