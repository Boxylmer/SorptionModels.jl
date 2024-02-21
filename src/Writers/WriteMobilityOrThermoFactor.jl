
DEFAULT_MOBILITY_THERMO_FACTOR_ANALYSIS_NAME = "Mobility and Thermo Factor Analysis"

"""
    write_analysis(analysis::MobilityFactorAnalysis, workbook::XLSX.XLSXFile; [name])
    write_analysis(analysis::MobilityFactorAnalysis, filepath::AbstractString; [name])
    write_analysis(analysis::ThermodynamicFactorAnalysis, workbook::XLSX.XLSXFile; [name])
    write_analysis(analysis::ThermodynamicFactorAnalysis, filepath::AbstractString; [name])

Write a ThermoFactorAnalysis or MobilityFactorAnalysis (which contains a thermodynamic factor analysis) out to a .xlsx `workbook`` or an .xlsx `filepath`. If a worksheet `name` is not provided, default names are created.

See [`ThermodynamicFactorAnalysis`](@ref) and [`MobilityFactorAnalysis`](@ref)
"""
function write_analysis(analysis::ThermodynamicFactorAnalysis, workbook::XLSX.XLSXFile; name=DEFAULT_MOBILITY_THERMO_FACTOR_ANALYSIS_NAME)
    if name in XLSX.sheetnames(workbook)
        sheet = workbook[name]
    else
        sheet = XLSX.addsheet!(workbook, name)
    end
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
    
    sheet[wrp, wcp + 1]  = "activity"
    sheet[wrp, wcp + 2]  = "activity σ"
    sheet[wrp, wcp + 3]   = "pressures"
    sheet[wrp, wcp + 4]   = "pressures σ"
    sheet[wrp, wcp + 5]   = "concentrations (cc/cc)"
    sheet[wrp, wcp + 6]   = "concentrations σ"
    sheet[wrp, wcp + 7]  = "ln(a)"
    sheet[wrp, wcp + 8]  = "ln(a) σ"
    sheet[wrp, wcp + 9]  = "mass fraction (ω)"
    sheet[wrp, wcp + 10]  = "ω σ"
    sheet[wrp, wcp + 11]  = "ln(ω)"
    sheet[wrp, wcp + 12]  = "ln(ω) σ"
    sheet[wrp, wcp + 13]  = "Thermodynamic Factor (α)"
    sheet[wrp, wcp + 14] = "α σ"
    
    wrp += 1
    
    sheet[wrp, wcp, dim=1] = collect(1:length(analysis.thermodynamic_factors))
    
    write_vector_of_maybe_measurements(wrp, wcp + 1, exp.(analysis.lna))
    write_vector_of_maybe_measurements(wrp, wcp + 3, analysis.pressures)
    write_vector_of_maybe_measurements(wrp, wcp + 5, analysis.concentrations)
    write_vector_of_maybe_measurements(wrp, wcp + 7, analysis.lna)
    write_vector_of_maybe_measurements(wrp, wcp + 9, exp.(analysis.lnw))
    write_vector_of_maybe_measurements(wrp, wcp + 11, analysis.lnw)
    write_vector_of_maybe_measurements(wrp, wcp + 13, analysis.thermodynamic_factors)
    return 15 # replace with whatever the next wcp should be (the mobility factor analysis will take off from there)
end

function write_analysis(analysis::MobilityFactorAnalysis, workbook::XLSX.XLSXFile; name=DEFAULT_MOBILITY_THERMO_FACTOR_ANALYSIS_NAME)
    sheet = XLSX.addsheet!(workbook, name)
    wrp = 1  # writer row position
    wcp = write_analysis(analysis.thermo_factor_analysis, workbook; name)  # writer col position

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
    
    sheet[wrp, wcp + 1]  = "kinetic factors (L) (cm2/s)"
    sheet[wrp, wcp + 2]  = "L σ"
    sheet[wrp, wcp + 3]  = "Diffusivity (D) (cm2/s)"
    sheet[wrp, wcp + 4]  = "D σ"

    wrp += 1
    
    write_vector_of_maybe_measurements(wrp, wcp + 1, analysis.kinetic_factors)
    write_vector_of_maybe_measurements(wrp, wcp + 3, analysis.kinetic_factors .* analysis.thermo_factor_analysis.thermodynamic_factors)
    
end

function write_analysis(analysis::Union{ThermodynamicFactorAnalysis, MobilityFactorAnalysis}, filepath::AbstractString; name=DEFAULT_MOBILITY_THERMO_FACTOR_ANALYSIS_NAME)
    if !isfile(filepath)
        XLSX.openxlsx(filepath, mode="w") do _ end
    end
    
    XLSX.openxlsx(filepath, mode="rw") do xf 
        return write_analysis(analysis, xf; name=name)
    end
end
