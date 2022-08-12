struct MobilityFactorAnalysis
    # will use this to hold resulting parameters from the deconvolution
    lna
    lnw
    slopes
    kinetic_factors
    thermodynamic_factors
end

"""
    MobilityFactorAnalysis(isotherm::IsothermData, transient_sorption_models::AbstractVector{<:TransientSorptionModel}, semi_thickness_cm::Number)

Deconvolute a vector of fitted transient sorption models, the corresponding equilibrium isotherm, and the semi-thickness of the polymer sample into a MobilityFactorAnalysis object.

# Arguments
- `isotherm::IsothermData`: Should be a single-component isotherm. If multiple components are present, only the first component will be deconvoluted.
* `transient_sorption_models::AbstractVector{<:TransientSorptionModel}`: Vector of TransientSorptionModel objects. Each fitted model needs to correspond to each step in the isotherm.
* `semi_thickness_cm::Number`: Half (semi) thickness of the polymer sample used in the sorption experiment in ``cm``.  
"""
function MobilityFactorAnalysis(isotherm::IsothermData, transient_sorption_models::AbstractVector{<:TransientSorptionModel}, semi_thickness_cm::Number)
    diffusivity_vector = get_diffusivity.(transient_sorption_models, semi_thickness_cm)
    return MobilityFactorAnalysis(isotherm, diffusivity_vector)
end

"""
    (isotherm::IsothermData, diffusivities::AbstractVector{<:Number})

Deconvolute an isotherm and already-known diffusivity values into their kinetic and thermodynamic components.

# Arguments
- `isotherm::IsothermData`: Should be a single-component isotherm. If multiple components are present, only the first component will be deconvoluted.
* `diffusivities::AbstractVector{<:Number}`: Vector of diffusivity values in ``cm^2/s``  
"""
function MobilityFactorAnalysis(isotherm::IsothermData, diffusivities::AbstractVector{<:Number})
    tfa = ThermodynamicFactorAnalysis(isotherm)
    kinetic_factors = tfa.slopes[1:length(diffusivities)].*diffusivities 
    return MobilityFactorAnalysis(tfa.lna, tfa.lnw, tfa.slopes, kinetic_factors, tfa.thermodynamic_factors)
end