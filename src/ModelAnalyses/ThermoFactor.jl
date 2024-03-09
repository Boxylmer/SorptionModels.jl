struct ThermodynamicFactorAnalysis
    # will use this to hold resulting parameters from the deconvolution
    pressures::Union{AbstractVector, Nothing}
    concentrations::AbstractVector
    lna::AbstractVector
    lnw::AbstractVector
    thermodynamic_factors::AbstractVector
end

"""
    ThermodynamicFactorAnalysis(isotherm::IsothermData)
Extract the thermodynamic factors of an isotherm using finite differencing approximations (no equations of state are used).

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
    penetrant_mass_fracs = penetrant_mass_fractions(isotherm; component=1)
    penetrant_activities = activities(isotherm, component=1)
    lna = log.(penetrant_activities)
    lnw = log.(penetrant_mass_fracs)
    slopes = estimate_slope_by_adjacent_points((lna), (lnw))
    thermo_factors = 1 ./ slopes
    pressures = partial_pressures(isotherm; component=1)
    if isnothing(pressures)
        fill!(similar(thermo_factors,String),"N/A")
    end

    concs = concentration(isotherm; component = 1)
    return ThermodynamicFactorAnalysis(pressures, concs, lna, lnw, thermo_factors)
end

"""
    ThermodynamicFactorAnalysis(isotherm, sorptionmodel, activity_function)
Extract the thermodynamic factors of an isotherm using analytical derivatives and exact activities from an equation of state (activity_function).

# Arguments
- `isotherm::IsothermData`: Should be a single-component isotherm. If multiple components are present, only the first component will be used.
- `sorptionmodel::SorptionModel`: SorptionModel that supports `predict_pressure`.
- `activity_function::Function`: Function which takes a pressure in MPa and returns an activity. 

Isotherms in this function will need pressures, polymer density, and penetrant molecular weight.
- Concentrations used in this method will be *predicted* concentrations from the `SorptionModel`, not tabulated concentrations. 
- The pressures used here will be found in the isotherm, so if you're using this with `MobilityFactorAnalysis`(@ref), you likely want to specify pressures manually. 

"""
function ThermodynamicFactorAnalysis(
    isotherm::IsothermData, 
    sorptionmodel::SorptionModel, 
    activity_function::Function, 
    )
    if isnothing(polymer_density(isotherm))
        throw(MissingException("Isotherm did not have polymer density."))
    end
    if isnothing(molecular_weights(isotherm))
        throw(MissingException("Isotherm did not have any molecular weights."))
    end
    ρ_pol = polymer_density(isotherm)
    mw_pen = molecular_weights(isotherm)[1]
    pressures = partial_pressures(isotherm; component=1)
    
    return ThermodynamicFactorAnalysis(pressures, ρ_pol, mw_pen, sorptionmodel, activity_function)
end

"""
    ThermodynamicFactorAnalysis(pressures, ρ_pol, mw_pen, sorptionmodel, activity_function)
Extract the thermodynamic factors of sorption data (see `ThermodynamicFactorAnalysis` method using an `isotherm`, `sorptionmodel`, `activity_function`).
This function allows for manual control of `pressures`, which is often necessary when using tabulated diffusivities rather than the isotherm's pressures. 

# Arguments
- `pressures::AbstractVector`: List of pressures (MPa) to evaluate the thermodynamic factor at.
- `ρ_pol::Number`: Polymer density in g/cm3.
- `mw_pwn::Number`: Penetrant molecular weight in g/mol. 
- `sorptionmodel::SorptionModel`: SorptionModel that supports `predict_pressure`.
- `activity_function::Function`: Function which takes a pressure in MPa and returns an activity. 
"""
function ThermodynamicFactorAnalysis(
        pressures::AbstractVector,
        ρ_pol::Number,
        mw_pen::Number,
        sorptionmodel::SorptionModel, 
        activity_function::Function, 
    )

    # # looking for ∂ln(a)/∂ln(ω)
    ln_a(pressure_mpa) = log(activity_function(pressure_mpa))
    ln_w(pressure_mpa) = log(
        ccpen_per_ccpol_to_mass_fractions(
            predict_concentration(sorptionmodel, pressure_mpa), ρ_pol, [mw_pen]
        )[2]
    )

    dlna_dp = ForwardDiff.derivative.(ln_a, pressures)
    dlnw_dp = ForwardDiff.derivative.(ln_w, pressures)

    α = dlna_dp ./ dlnw_dp

    lna = ln_a.(pressures)
    lnw = ln_w.(pressures)

    
    concs =  predict_concentration(sorptionmodel, pressures)

    return ThermodynamicFactorAnalysis(pressures, concs, lna, lnw, α)
end
