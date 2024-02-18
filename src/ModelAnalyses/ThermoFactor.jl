struct ThermodynamicFactorAnalysis
    # will use this to hold resulting parameters from the deconvolution
    lna
    lnw
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
    return ThermodynamicFactorAnalysis(lna, lnw, thermo_factors)
end


"""
    ThermodynamicFactorAnalysis(isotherm::IsothermData)
Extract the thermodynamic factors of an isotherm's steps.

# Arguments
- `isotherm::IsothermData`: Should be a single-component isotherm. If multiple components are present, only the first component will be used.
- `sorptionmodel::SorptionModel`: SorptionModel that supports `predict_pressure`.
- `activity_function::Function`: Function which takes a pressure in MPa and returns an activity. 

Isotherms in this function will need polymer density, penetrant molecular weight.
"""
function ThermodynamicFactorAnalysis(
        isotherm::IsothermData, 
        sorptionmodel::SorptionModel, 
        activity_function::Function, 
    )

    if isnothing(polymer_density(isotherm))
        throw(MissingException("Isotherm did not have polymer density."))
    end
    if isnothing(penetrant_molecular_weight(isotherm))
        throw(MissingException("Isotherm did not have any molecular weights."))
    end

    # looking for ∂ln(a)/∂ln(ω)
    function ln_a(ln_w)
        mass_faction = exp(ln_w)
        conc = polymer_phase_mass_fractions_to_ccpen_per_ccpol([1-mass_faction, mass_faction], polymer_density(isotherm), [penetrant_molecular_weight])
        pres = predict_pressure(sorptionmodel, conc)
        activity = activity_function(pres)
        return ln(activity)
    end

    ln_a(pressure_mpa) = log(activity_function(pressure_mpa))
    ln_w(pressure_mpa) = log(
        ccpen_per_ccpol_to_mass_fractions(
            predict_concentration(sorptionmodel, pressure_mpa), polymer_density, [penetrant_molecular_weight]
        )
    )

    dlna_dp = ForwardDiff.derivative.(ln_a, partial_pressures(isotherm; component=1))
    dlnw_dp = ForwardDiff.derivative.(ln_w, partial_pressures(isotherm; component=1))

    α = dlna_dp ./ dlnw_dp

    lna = ln_a(partial_pressures(isotherm; component=1))
    lnw = ln_w(partial_pressures(isotherm; component=1))

    # single component
    return ThermodynamicFactorAnalysis(lna, lnw, α)
end
