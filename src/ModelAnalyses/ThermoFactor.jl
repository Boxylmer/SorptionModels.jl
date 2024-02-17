struct ThermodynamicFactorAnalysis
    # will use this to hold resulting parameters from the deconvolution
    lna
    lnw
    slopes
    thermodynamic_factors
end

"""
    ThermodynamicFactorAnalysis(isotherm::IsothermData)
Extract the thermodynamic factors of an isotherm's steps.

# Arguments
- `isotherm::IsothermData`: Should be a single-component isotherm. If multiple components are present, only the first component will be used.

Isotherms in this function will need penetrant mass fractions and activities.
"""

function ThermodynamicFactorAnalysis(isotherm::IsothermData)
    if isnothing(activities(isotherm))
        throw(MissingException("Isotherm did not have any activities."))
    end
    if isnothing(penetrant_mass_fractions(isotherm))
        throw(MissingException("Isotherm did not have any mass fractions."))
    end

    # single component
    penetrant_mass_fracs = penetrant_mass_fractions(isotherm; component=1)  # component 1 
    penetrant_activities = activities(isotherm, component=1)
    lna = log.(penetrant_activities)
    lnw = log.(penetrant_mass_fracs)
    slopes = estimate_slope_by_adjacent_points((lna), (lnw))
    thermo_factors = 1 ./ slopes
    return ThermodynamicFactorAnalysis(lna, lnw, slopes, thermo_factors)
end