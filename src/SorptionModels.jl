module SorptionModels

    using MembraneBase
    using MembraneEOS
    using Optim
    using Measurements

    include(joinpath("SorptionModels", "SorptionModels.jl"))
    export predict_concentration
    export fit_model

    include(joinpath("SorptionModels", "DualMode.jl"))
    export fit_dualmode_model
    
    include(joinpath("SorptionModels", "GAB.jl"))
    export fit_gab_model

    include(joinpath("SorptionModels", "NELF.jl"))
    export NELFModel
    export calculate_polymer_phase_density
    export calculate_swelled_polymer_density


end
