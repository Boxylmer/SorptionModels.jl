using SorptionModels
using Measurements
using MembraneBase
using MembraneEOS
using SorptionApparatus
using XLSX 

function ZimmLundberg_Gas_Analysis(sample_path::String,output_path::String, pen_vapor_pressure, pen_molar_volume; sorption_model = :DualMode)
    sample_xlsx = sample_path
    sample_output_xlsx = output_path

    sample_sorption_model = processtemplate(
        GasSorptionApparatus(), 
        sample_xlsx, 
        sample_output_xlsx; 
        overwrite=true)


    function pressure_to_act(pressure)
        act = pressure/pen_vapor_pressure
        return act
    end
    
    isotherm_pressures = partial_pressures(sample_sorption_model.isotherm; component = 1)
    activities = pressure_to_act.(isotherm_pressures)
    concs = concentration(sample_sorption_model.isotherm; component = 1)


    if sorption_model == :DualMode
        sample_dm = fit_model(
            DualMode(),
            sample_sorption_model.isotherm; 
            use_fugacity = true, 
            uncertainty_method = :JackKnife)
       
        zimm_lundberg_analysis = ZimmLundbergAnalysis(
            sample_dm,
            isotherm_pressures, 
            pen_vapor_pressure, 
            pen_molar_volume)
    elseif sorption_model == :GAB  
        sample_GAB = fit_model(
            GAB(),
            activities,
            concs;
            uncertainty_method = :JackKnife)
        
        zimm_lundberg_analysis = ZimmLundbergAnalysis(
            sample_GAB,
            activities,
            pen_vapor_pressure,
            pen_molar_volume)
    end 

    # @show sample_GAB.cp
    # @show sample_GAB.k
    # @show sample_GAB.a
    
    XLSX.openxlsx(sample_output_xlsx, mode = "rw") do xf
        write_analysis(zimm_lundberg_analysis, xf)
    end
end



sample_path = joinpath(@__DIR__,"CPIM_NH3_sorption.xlsx")
output_path_GAB = joinpath(@__DIR__,"diagnostics_output","CPIM_NH3_sorption_output_GAB.xlsx")
output_path_DM = joinpath(@__DIR__,"diagnostics_output","CPIM_NH3_sorption_output_DM.xlsx")

ammonia_molar_volume_35C = 28.956 #cm^3/mol
ammonia_vapor_pressure_35C = 1.35 #MPa 

function ammonia_pressure_to_activity(pressure::Number)
    pvap = 1350/1000 #MPa
    activity = pressure / pvap
    return activity
end

ZimmLundberg_Gas_Analysis(sample_path,output_path_GAB, ammonia_vapor_pressure_35C, ammonia_molar_volume_35C; sorption_model = :GAB)
ZimmLundberg_Gas_Analysis(sample_path,output_path_DM, ammonia_vapor_pressure_35C, ammonia_molar_volume_35C)
#



