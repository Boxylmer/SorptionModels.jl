module SorptionModels

    using MembraneBase
    using MembraneEOS
    using ForwardDiff
    using Optim
    using Measurements
    using StaticArrays

    using XLSX
    using DelimitedFiles

    using Plots # for diagnostic methods todo remove this


    include("ForwardDiffOverMeasurements.jl")

    include(joinpath("SorptionModels", "ModelMethods.jl"))
    export predict_concentration
    export a_predict_concentration
    export predict_pressure, predict_activity
    export infinite_dilution_solubility
    export fit_model

    include(joinpath("SorptionModels", "DualMode.jl"))
    export fit_dualmode_model
    export DualMode, DualModeModel
    
    include(joinpath("SorptionModels", "GAB.jl"))
    export fit_gab_model
    export GAB, GABModel

    include(joinpath("SorptionModels", "NELF.jl"))
    export NELF, NELFModel
    export calculate_polymer_phase_density
    export calculate_swelled_polymer_density
    export fit_kij, fit_ksw

    include(joinpath("SorptionModels", "DGRPT.jl"))
    export DGRPT, DGRPTModel

    include(joinpath("TransientSorptionModels", "TransientSorptionModels.jl"))
    export fit_transient_sorption_model
    export get_diffusivity
    export predict_sorption

    include(joinpath("TransientSorptionModels", "FickianSorption.jl"))
    export FickianSorption, FickianSorptionModel
    include(joinpath("TransientSorptionModels", "BerensHopfenbergSorption.jl"))
    export BerensHopfenbergSorption, BerensHopfenbergSorptionModel
    include(joinpath("TransientSorptionModels", "ModifiedBerensHopfenbergSorption.jl"))
    export ModifiedBerensHopfenbergSorption, ModifiedBerensHopfenbergSorptionModel


    # Model specific analysis
    include(joinpath("ModelAnalyses", "VantHoffDualMode.jl"))
    export VantHoffDualModeAnalysis

    include(joinpath("ModelAnalyses", "IsostericHeatOfSorption.jl"))
    export IsostericHeatAnalysis

    include(joinpath("ModelAnalyses", "ThermoFactor.jl"))
    export ThermodynamicFactorAnalysis

    include(joinpath("ModelAnalyses", "MobilityFactor.jl"))
    export MobilityFactorAnalysis

    include(joinpath("ModelAnalyses", "PartialImmobilization.jl"))
    export PartialImmobilizationModel

    include(joinpath("ModelAnalyses", "ZimmLundberg.jl"))
    export ZimmLundbergAnalysis

    include(joinpath("ModelAnalyses", "DualModeDesorption.jl"))
    export DualModeDesorption

    # model writers
    include(joinpath("Writers", "WriteVantHoffDualMode.jl"))
    include(joinpath("Writers", "WriteIsostericHeat.jl"))
    # mobility factor not implemented yet
    # thermodynamic factor not implemented yet
    include(joinpath("Writers", "WritePartialImmobilization.jl"))
    include(joinpath("Writers", "WriteZimmLundberg.jl"))
    include(joinpath("Writers", "WriteDualModeDesorption.jl"))
    export write_analysis


    # Diagnostic Methods
    include(joinpath("ModelDiagnostics", "NELFDiagnostics.jl"))
    export nelf_characteristic_parameter_error_map
end
