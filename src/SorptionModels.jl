module SorptionModels

    using MembraneBase
    using MembraneEOS
    using Optim
    using Measurements

    export predict_concentration
    export fit_model

    include(joinpath("SorptionModels", "SorptionModels.jl"))

    include(joinpath("SorptionModels", "DualMode.jl"))
    export fit_dualmode_model
    
    include(joinpath("SorptionModels", "GAB.jl"))
    export fit_gab_model

    include(joinpath("SorptionModels", "NELF.jl"))
    export NELFModel
    export calculate_polymer_phase_density
    export calculate_swelled_polymer_density

    include(joinpath("TransientSorptionModels", "TransientSorptionModels.jl"))
    export fit_transient_sorption_model
    export get_diffusivity
    export predict_sorption

    include(joinpath("TransientSorptionModels", "FickianSorption.jl"))
    export FickianSorptionModel
    include(joinpath("TransientSorptionModels", "BerensHopfenbergSorption.jl"))
    export BerensHopfenbergSorptionModel
    include(joinpath("TransientSorptionModels", "ModifiedBerensHopfenbergSorption.jl"))
    export ModifiedBerensHopfenbergSorptionModel


end
