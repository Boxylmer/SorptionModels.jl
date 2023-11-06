
module IsostericHeatAnalysisHelper
    const default_isosteric_heat_sheet_name = "Isosteric Heat Analysis"
end
IHAH = IsostericHeatAnalysisHelper

"""
    write_analysis(analysis::AbstractIsostericHeatAnalysis, workbook::XLSX.XLSXFile; [name])
    write_analysis(analysis::AbstractIsostericHeatAnalysis, filepath::AbstractString; [name])

Write an isosteric heat `analysis` out to a .xlsx `workbook`` or an .xlsx `filepath`. If a worksheet `name` is not provided, default names are created.
See [`IsostericHeatAnalysis`](@ref)
"""
function write_analysis(analysis::AbstractIsostericHeatAnalysis, workbook::XLSX.XLSXFile; name=IHAH.default_isosteric_heat_sheet_name)
    sheet = XLSX.addsheet!(workbook, name)

    num_isotherms = length(analysis.isotherms)
    writer_row_position = 1
    writer_col_position = 1

    # first show the fitted, interpolated data
    sheet[writer_row_position, writer_col_position] = "Fitted, interpolated data" 
    writer_row_position += 1
    sheet[writer_row_position, writer_col_position] = "Sampled CC/CC"
    sheet[writer_row_position, writer_col_position + 1, dim=2] = [measurement(conc).val for conc in analysis.sampled_concentrations]
    writer_row_position +=1
    sheet[writer_row_position, writer_col_position] = "T"
    sheet[writer_row_position + 1, writer_col_position, dim=1] = [measurement(temp).val for temp in analysis.temperatures]
    for (idx, pressure_vector) in enumerate(analysis.pressure_vectors)
        sheet[writer_row_position, writer_col_position + idx] = "P[MPa]"
        sheet[writer_row_position + 1, writer_col_position + idx, dim=1] = [measurement(meas).val for meas in pressure_vector]
    end

    # next show the resulting ln(P) vs 1/T curves 
    writer_row_position += num_isotherms + 2

    if typeof(analysis) <: IsostericHeatAnalysis
        regression_units = "Manipulated data for getting ∂ln(P)/∂(1/T)" 
    elseif typeof(analysis) <: WebbIsostericHeatAnalysis
        regression_units = "Manipulated data for getting ∂ln(L(STP)/mol / V[L/mol])/∂(1/T)" # TODO check units for 22.313 L mol
    end
    sheet[writer_row_position, writer_col_position] = regression_units
    writer_row_position += 1
    sheet[writer_row_position, writer_col_position] = "Sampled CC/CC"
    sheet[writer_row_position, writer_col_position + 1, dim=2] = [measurement(conc).val for conc in analysis.sampled_concentrations]
    writer_row_position +=1
    sheet[writer_row_position, writer_col_position] = "1/T (K^-1)"
    sheet[writer_row_position + 1, writer_col_position, dim=1] = [measurement(temp).val for temp in analysis.inverse_temperature_vector]
    
    if typeof(analysis) <: IsostericHeatAnalysis
        for (idx, vec) in enumerate(analysis.ln_pressure_vectors)
            sheet[writer_row_position, writer_col_position + idx] = "ln(P[MPa])"
            sheet[writer_row_position + 1, writer_col_position + idx, dim=1] = [measurement.(vec_meas).val for vec_meas in vec] 
        end
    elseif typeof(analysis) <: WebbIsostericHeatAnalysis
        for (idx, vec) in enumerate(analysis.ln_inv_vol_curves)
            sheet[writer_row_position, writer_col_position + idx] = "ln(L(STP)/mol / V[L/mol])" 
            sheet[writer_row_position + 1, writer_col_position + idx, dim=1] = [measurement.(vec_meas).val for vec_meas in vec] 
        end
    end
   

    # finally, show the found isosteric heats as a function of concentration
    writer_row_position += num_isotherms + 2
    sheet[writer_row_position, writer_col_position] = "Isosteric heats as a function of concentration" 
    writer_row_position += 1
    sheet[writer_row_position, writer_col_position] = "Conc (CC/CC)" 
    sheet[writer_row_position, writer_col_position + 1] = "ΔH_s (J/mol)" 
    sheet[writer_row_position, writer_col_position + 2] = "ΔH_s err (J/mol)" 
    sheet[writer_row_position, writer_col_position + 3] = "ΔS_s (J/molK)" 
    sheet[writer_row_position, writer_col_position + 4] = "ΔS_s err (J/molK)" 
    sheet[writer_row_position, writer_col_position + 5] = "S0 (J/molK) (S0 = exp(ΔS_s/R))" 
    sheet[writer_row_position, writer_col_position + 6] = "S0 err (J/molK)"

    writer_row_position += 1
    sheet[writer_row_position, writer_col_position, dim=1] = [measurement(conc).val for conc in analysis.sampled_concentrations]
    sheet[writer_row_position, writer_col_position + 1, dim=1] = [measurement(heat).val for heat in analysis.isosteric_heat_at_conc]
    sheet[writer_row_position, writer_col_position + 2, dim=1] = [measurement(heat).err for heat in analysis.isosteric_heat_at_conc]

    sheet[writer_row_position, writer_col_position + 3, dim=1] = [measurement(entropy).val for entropy in analysis.isosteric_entropy_at_conc]
    sheet[writer_row_position, writer_col_position + 4, dim=1] = [measurement(entropy).err for entropy in analysis.isosteric_entropy_at_conc]

    sheet[writer_row_position, writer_col_position + 5, dim=1] = [measurement(factor).val for factor in analysis.pre_exponential_factors]
    sheet[writer_row_position, writer_col_position + 6, dim=1] = [measurement(factor).err for factor in analysis.pre_exponential_factors]
end

# todo all of these "write analysis" functions, when given a filepath, could be combined
function write_analysis(analysis::AbstractIsostericHeatAnalysis, filepath::AbstractString; name=missing)
    if ismissing(name)
        if typeof(analysis) <: IsostericHeatAnalysis
            name = IHAH.default_isosteric_heat_sheet_name
        elseif typeof(analysis) <: WebbIsostericHeatAnalysis
            name = "WebbIsostericHeatAnalysis"
        else
            name = "Abstract or Custom Isosteric Heat Analysis"
        end
    end
    
    if !isfile(filepath)
        XLSX.openxlsx(filepath, mode="w") do _ end
    end

    XLSX.openxlsx(filepath, mode="rw") do xf  
        return write_analysis(analysis, xf; name=name)
    end
end
