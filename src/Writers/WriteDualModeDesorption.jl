
module DualModeDesorptionHelper
    const default_sheet_name = "Dual Mode Desorption Analysis"
end
DMDH = DualModeDesorptionHelper




"""
    write_analysis(analysis::DualModeDesorption, workbook::XLSX.XLSXFile; [name])
    write_analysis(analysis::DualModeDesorption, filepath::AbstractString; [name])

Write a Dual Mode Desorption `analysis` out to a .xlsx `workbook`` or an .xlsx `filepath`. If a worksheet `name` is not provided, default names are created.

See [`DualModeDesorption`](@ref)
"""
function write_analysis(analysis::DualModeDesorption, workbook::XLSX.XLSXFile; name=DMDH.default_sheet_name, npts=30)
    sheet = XLSX.addsheet!(workbook, name)

    wrp = 1  # writer row position
    wcp = 1 # writer column position

    if using_fugacity(analysis)
        p_text = "fugacity (MPa)"
        original_ps = strip_measurement_to_value(fugacity(analysis.isotherm; component=1))
    else
        p_text = "pressure (MPa)"
        original_ps = strip_measurement_to_value(partial_pressures(analysis.isotherm; component=1))
    end
    original_cs = concentration(analysis.isotherm; component=1)

    ps = collect(range(0, analysis.max_p, npts))
    sorption_cs = predict_concentration(analysis.sorbing_model, ps)
    desorption_cs = predict_concentration(analysis.desorbing_model, ps)

    function write_model_to_row(m::DualModeModel, rowstart, colstart)
        sheet[rowstart, colstart] = measurement(m.ch).val
        sheet[rowstart, colstart+1] = measurement(m.ch).err
        sheet[rowstart, colstart+2] = measurement(m.b).val
        sheet[rowstart, colstart+3] = measurement(m.b).err
        sheet[rowstart, colstart+4] = measurement(m.kd).val
        sheet[rowstart, colstart+5] = measurement(m.kd).err
    end

    # header
    sheet[wrp, wcp, dim=2] = ["Fitted to $p_text", "CH' (CC/CC)", "err", "b (1/MPa)", "err", "kd ((CC/CC)/MPa)", "err"]

    wrp += 1
    sheet[wrp, wcp, dim=1] = ["Sorption", "Desorption", "% Change"]

    write_model_to_row(analysis.sorbing_model, wrp, wcp + 1)
    wrp += 1
    write_model_to_row(analysis.desorbing_model, wrp, wcp + 1)

    wrp += 1
    wcp += 1

    percent_change_ch = (analysis.sorbing_model.ch - analysis.desorbing_model.ch) / analysis.sorbing_model.ch * 100
    percent_change_b = (analysis.sorbing_model.b - analysis.desorbing_model.b) / analysis.sorbing_model.b * 100
    percent_change_kd = (analysis.sorbing_model.kd - analysis.desorbing_model.kd) / analysis.sorbing_model.kd * 100
    sheet[wrp, wcp] = measurement(percent_change_ch).val
    sheet[wrp, wcp+1] = measurement(percent_change_ch).err
    sheet[wrp, wcp+2] = measurement(percent_change_b).val
    sheet[wrp, wcp+3] = measurement(percent_change_b).err
    sheet[wrp, wcp+4] = measurement(percent_change_kd).val
    sheet[wrp, wcp+5] = measurement(percent_change_kd).err

    wrp += 3
    wcp = 1
    
    sheet[wrp, wcp, dim=2] = [p_text, "Sorption CC/CC", "err", "Desorption CC/CC", "err", "Original $p_text", "Original concentration (CC/CC)", "err"] 
    wrp += 1
    sheet[wrp, wcp, dim=1] = ps
    wcp += 1

    function write_column_maybe_uncertainty(wrp, wcp, vec)
        for (idx, maybe_meas) in enumerate(vec) 
            sheet[wrp + idx - 1, wcp] = measurement(maybe_meas).val
            sheet[wrp + idx - 1, wcp + 1] = measurement(maybe_meas).err
        end
    end

    write_column_maybe_uncertainty(wrp, wcp, sorption_cs)
    wcp += 2
    
    write_column_maybe_uncertainty(wrp, wcp, desorption_cs)
    wcp += 2

    sheet[wrp, wcp, dim=1] = original_ps
    wcp += 1
    write_column_maybe_uncertainty(wrp, wcp, original_cs)

end

# todo all of these "write analysis" functions, when given a filepath, could be combined
function write_analysis(analysis::DualModeDesorption, filepath::AbstractString; name=DMDH.default_sheet_name)
    if !isfile(filepath)
        XLSX.openxlsx(filepath, mode="w") do _ end
    end

    XLSX.openxlsx(filepath, mode="rw") do xf  
        return write_analysis(analysis, xf; name=name)
    end
end
