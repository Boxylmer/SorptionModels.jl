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

Deconvolute a vector of fitted transient sorption models, the corresponding equilibrium isotherm, and the semi-thickness of the polymer sample into a DiffusivityDeconvolution object.

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
    if num_components(isotherm) > 1
        throw(ErrorException("Deconvoluting diffusivity isn't implemented for multicomponent isotherms yet"))
    end
    if !(num_steps(isotherm) == length(diffusivities))
        throw(DimensionMismatch("Lengths of isotherm steps and diffusion coefficients didn't line up."))
    end
    if isnothing(activities(isotherm))
        throw(MissingException("Isotherm did not have any activities"))
    end
    # if isnothing(molecular_weights(isotherm))
    #     throw(MissingException("Isotherm does not contain molecular weights"))
    # end

    # single component
    penetrant_mass_fracs = penetrant_mass_fractions(isotherm; component=1)  # component 1 
    penetrant_activities = activities(isotherm, component=1)
    lna = log.(penetrant_activities)
    lnw = log.(penetrant_mass_fracs)
    slopes = estimate_slope_by_adjacent_points.(Ref(lna), Ref(lnw), 1:length(lna))

    kinetic_factors = slopes[1:length(diffusivities)].*diffusivities 
    thermo_factors = 1 ./ slopes
    
    return MobilityFactorAnalysis(lna, lnw, slopes, kinetic_factors, thermo_factors)
end