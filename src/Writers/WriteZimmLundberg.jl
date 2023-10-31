"""
    write_analysis(analysis::ZimmLundbergAnalysis, workbook::XLSX.XLSXFile; [name])
    write_analysis(analysis::ZimmLundbergAnalysis, filepath::AbstractString; [name])

Write a Zimm Lundberg `analysis` out to a .xlsx `workbook`` or an .xlsx `filepath`. If a worksheet `name` is not provided, default names are created.

See [`ZimmLundbergAnalysis`](@ref)
"""
function write_analysis(analysis::ZimmLundbergAnalysis, workbook::XLSX.XLSXFile; name="Zimm Lundberg Analysis")
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
        sheet[row_start + rowpos, col_start + colpos] = "Step"
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

    sheet[wrp, wcp] = "Step"
    sheet[wrp, wcp + 1] = "Activity"
    sheet[wrp, wcp + 2] = "Activity σ"
    sheet[wrp, wcp + 3] = "Volume Fraction"
    sheet[wrp, wcp + 4] = "Volume Fraction σ"
    sheet[wrp, wcp + 5] = "a/φ"
    sheet[wrp, wcp + 6] = "a/φ σ"
    sheet[wrp, wcp + 7] = "∂(a/φ)/∂(a)"
    sheet[wrp, wcp + 8] = "∂(a/φ)/∂(a) σ"
    sheet[wrp, wcp + 9] = "Cluster Function"
    sheet[wrp, wcp + 10] = "Cluster Function σ"
    sheet[wrp, wcp + 11] = "Avg. Cluster Size"
    sheet[wrp, wcp + 12] = "Avg. Cluster Size σ"
    
    wrp+=1
    
    nsteps = length(analysis.activities)
    write_vector_of_maybe_measurements(wrp, wcp, 1:1:nsteps; write_errors=false)
    wcp+=1
    write_vector_of_maybe_measurements(wrp, wcp, analysis.activities)
    wcp+=2
    write_vector_of_maybe_measurements(wrp, wcp, analysis.vol_fracs) 
    wcp+=2
    write_vector_of_maybe_measurements(wrp, wcp, analysis.a_over_phi_values)
    wcp+=2
    write_vector_of_maybe_measurements(wrp, wcp, analysis.a_over_phi_derivatives)
    wcp+=2
    write_vector_of_maybe_measurements(wrp, wcp, analysis.cluster_functions)
    wcp+=2
    write_vector_of_maybe_measurements(wrp, wcp, analysis.average_cluster_size)

end

function write_analysis(analysis::ZimmLundbergAnalysis, filepath::AbstractString; name="Zimm Lundberg Analysis")
    if !isfile(filepath)
        XLSX.openxlsx(filepath, mode="w") do _ end
    end
    
    XLSX.openxlsx(filepath, mode="rw") do xf 
        return write_analysis(analysis, xf; name=name)
    end
end
