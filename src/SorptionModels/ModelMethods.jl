abstract type SorptionModel end

"""
    predict_concentration(model::SorptionModel, args...)

Predict concentration with a sorption model. Predictions are generally returned in the format of a vector of values corresponding to an input of a vector of pressures.
"""
function predict_concentration end

"""
    fit_model(model, ::IsothermData)

Fit a sorption model to an isotherm. If the data required for the model is not in the isotherm, an error message will be returned.

Any keyword arguments get passed on to the specific model.
"""
function fit_model(model::SorptionModel, ::IsothermData) end

#TODO: this could be useful in Clapeyron?
function VT_chemical_potential_at_i(model,V,T,z = Clapeyron.SA[1.0],i)
    f(_z) = Clapeyron.eos(model,V,T,_z)
    return Clapeyron.Solvers.grad_at_i(f,z,i)
end

function VT_chemical_potential_res_at_i(model,V,T,z = Clapeyron.SA[1.0],i)
    f(_z) = Clapeyron.eos_res(model,V,T,_z)
    return Clapeyron.Solvers.grad_at_i(f,z,i)
end

#chemical potential in terms of mass density (g/cm3), mass fractions
"""
    ρTw_chemical_potential(model,ρ,T,w)

returns the chemical potential of a mixture model.
inputs:
- `model`: `Clapeyron.EoSModel`
- `ρ`: density [g/cm3]
- `T`: temperature [K]
- `w`: vector of mass fractions [kg/kg]

returns:
 - `μ` vector of chemical potentials [J/mol]
"""
function ρTw_chemical_potential(model,ρ,T,w)
    #vector of molecular weights, g/mol
    MW = Clapeyron.mw(model)
    #mol fractions
    n = MembraneBase.mass_fractions_to_mole_fractions(w,MW)
    m = molar_mass(model,n) #kg/mol
    #1 g/cm3 = 1000 kg/m3
    ρ_SI = ρ*1000
    V = m/ρ_SI
    return Clapeyron.VT_chemical_potential(model,V,T,n)
end

"""
    ub_density(model,w,input = :molar)

returns the upper bound density, in [g/cm3].
it calls `Clapeyron.lb_volume`.
inputs:
- `model`: `Clapeyron.EoSModel`
- `w`: vector of mass fractions [kg/kg]
- `input` (optional): a symbol specifying if `x` is in molar or mass base


returns:
 - `ub_rho` upper bound density [g/cm3]
"""
function ub_density(model,x,input = :molar) 
    if input in (:molar,:mol,:n)
        MW = Clapeyron.mw(model)
        n = MembraneBase.mass_fractions_to_mole_fractions(x,MW)
        lb_v = Clapeyron.lb_volume(model,n) #m3/mol
    else
        lb_v = Clapeyron.lb_volume(model,x)
    end
    m = molar_mass(model,x,input)
    ub_rho_SI = m/lb_v
    ub_rho = ub_rho_SI/1000 #(1 kg/m3 = 0.001 g/cm3)
    return ub_rho
end

"""
    molar_mass(model,x,input::Symbol = :molar)

returns the amount of kg in a mol of mixture [kg/mol]
if `input == :molar`, then `x` is assumed to be molar fractions. if `input == :weight`, then `x` is assumed to be mass fractions.
inputs:
- `model`: `Clapeyron.EoSModel`
- `x`: vector of molar or mass compositions
- `input` (optional): a symbol specifying if `x` is in molar or mass base

returns:
 - `m`: molar mass [kg/mol]
"""
function molar_mass(model::Clapeyron.EoSModel,x,input = :molar)
    if input in (:molar,:mol,:n)
        return Clapeyron.molecular_weight(model,x)
    else
        MW = Clapeyron.mw(model)
        m = 0.0*zero(eltype(x))
        for i in eachindex(x)
            wi,mwi = x[i],MW[i]
            #kgi/kgt / (gi/moli)
            #kgi/kgt / moli/gi * 1000/1 gi/kgi
            #moli/kgt
            m += 1000*wi/mwi
        end
        return 1/m
    end
end

"""
    ρTw_chemical_potential_at_i(model,ρ,T,w,i)

returns the chemical potential of a mixture model at the index i.
inputs:
- `model`: `Clapeyron.EoSModel`
- `ρ`: density [g/cm3]
- `T`: temperature [K]
- `w`: vector of mass fractions [kg/kg]
- `i`: index for requested chemical potential
returns:
 - `μi` chemical potential at index i [J/mol]
"""
function ρTw_chemical_potential_at_i(model,ρ,T,w,i)
    #vector of molecular weights, g/mol
    MW = Clapeyron.mw(model)
    #mol fractions
    n = MembraneBase.mass_fractions_to_mole_fractions(w,MW)
    m = molar_mass(model,n) #kg/mol
    #1 g/cm3 = 1000 kg/m3
    ρ_SI = ρ*1000
    V = m/ρ_SI
    return VT_chemical_potential_at_i(model,V,T,n,i)
end
