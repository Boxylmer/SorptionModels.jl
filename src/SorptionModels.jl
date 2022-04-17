module SorptionModels

    using MembraneBase
    using MembraneEOS
    using Optim
    using Measurements

    using XLSX
    using DelimitedFiles

    include(joinpath("SorptionModels", "ModelMethods.jl"))
    export predict_concentration
    export fit_model

    include(joinpath("SorptionModels", "DualMode.jl"))
    export fit_dualmode_model
    export DualMode, DualModeModel
    
    include(joinpath("SorptionModels", "GAB.jl"))
    export fit_gab_model
    export GAB, GABModel

    include(joinpath("SorptionModels", "NELF.jl"))
    export NELF
    export NELFModel
    export calculate_polymer_phase_density
    export calculate_swelled_polymer_density

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

    include(joinpath("ModelAnalyses", "MobilityFactor.jl"))
    export MobilityFactorAnalysis

    include(joinpath("ModelAnalyses", "PartialImmobilization.jl"))
    export PartialImmobilizationModel

    # model writers
    include(joinpath("Writers", "WriteVantHoffDualMode.jl"))
    include(joinpath("Writers", "WriteIsostericHeat.jl"))
    # mobility factor not implemented yet
    include(joinpath("Writers", "WritePartialImmobilization.jl"))
    export write_analysis

end
