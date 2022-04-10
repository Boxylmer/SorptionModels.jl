struct DGRPT end
const default_dgrpt_taylor_expansion_order = 12

"""

    Requires that in the kijmatrix, the polymer is the first index. The remaining indexes must match `penetrants`.

"""

struct DGRPTModel{BMT, POLYMT, PDT} <: SorptionModel
    bulk_model::BMT                 # [ChemicalParameters]
    polymer_model::POLYMT           # ChemicalParameters
    polymer_dry_density::PDT        # number
end

function predict_concentration(
        model::DGRPTModel, temperature, pressure, bulk_penetrant_mole_fractions, penetrant_molecular_weights=nothing; 
        taylor_series_order=default_dgrpt_taylor_expansion_order, units=:cc)
    penetrant_mass_fraction_initial_guesses = ones(length(bulk_penetrant_mole_fractions)) * 1e-5
   
    # Optim.jl target
    mass_fraction_error = make_penetrant_mass_fraction_target(model, temperature, pressure, bulk_penetrant_mole_fractions; taylor_series_order)

    lower = ones(length(penetrant_mass_fraction_initial_guesses))*eps()
    upper = ones(length(penetrant_mass_fraction_initial_guesses)) .- eps()
    res = Optim.optimize(
        mass_fraction_error, lower, upper, 
        penetrant_mass_fraction_initial_guesses, 
        Fminbox(LBFGS()),
        Optim.Options(
            allow_f_increases = false,
            x_tol = 1e-4,
            # g_tol = 1e-7,
            # f_tol = 1e-4,
        ); autodiff=:forward)
    penetrant_mass_fractions = Optim.minimizer(res)
    # end Optim.jl

    polymer_phase_mass_fractions = vcat(1 - sum(penetrant_mass_fractions), penetrant_mass_fractions)
    if units==:frac
        return polymer_phase_mass_fractions
    elseif units==:g
        concs_g_g = polymer_phase_mass_fractions_to_gpen_per_gpol(polymer_phase_mass_fractions)
        return concs_g_g
    elseif units==:cc
        if isnothing(penetrant_molecular_weights)
            throw(ArgumentError("Need to specify penetrant molecular weights (g/mol) to get concentration in units of cc/cc!"))
        end
        concs_cc_cc = polymer_phase_mass_fractions_to_ccpen_per_ccpol(polymer_phase_mass_fractions, model.polymer_dry_density, penetrant_molecular_weights)
        return concs_cc_cc
    end
end

function polymer_density(model::DGRPTModel, temperature, pressure, bulk_penetrant_mole_fractions, penetrant_molecular_weights; taylor_series_order=default_dgrpt_taylor_expansion_order)
    polymer_phase_mass_fractions = predict_concentration(model, temperature, pressure, bulk_penetrant_mole_fractions, penetrant_molecular_weights; taylor_series_order, units=:frac)
    polymer_density = solve_polymer_density(model, temperature, polymer_phase_mass_fractions; taylor_series_order)
    return polymer_density
end

function calculate_bulk_phase_chemical_potentials(model::DGRPTModel, temperature, pressure, bulk_phase_mole_fractions)
    μ = chemical_potential(
            model.bulk_model, 
            pressure, 
            temperature, 
            bulk_phase_mole_fractions)
    return μ
end

function calculate_bulk_phase_activities(model::DGRPTModel, temperature, pressure, bulk_phase_mole_fractions)
    a = activity(
            model.bulk_model, 
            pressure, 
            temperature, 
            bulk_phase_mole_fractions)
    return a
end

function calculate_polymer_phase_chemical_potentials(model::DGRPTModel, temperature, polymer_density, polymer_phase_mass_fractions)
    polymer_phase_density = polymer_density / polymer_phase_mass_fractions[1]
    μ = ρTω_chemical_potential(
        model.polymer_model, 
        polymer_phase_density, 
        temperature, 
        polymer_phase_mass_fractions)
    return μ
end

function calculate_polymer_phase_activities(model::DGRPTModel, temperature, polymer_density, polymer_phase_mass_fractions)
    polymer_phase_density = polymer_density / polymer_phase_mass_fractions[1]
    a = ρTω_activity(
        model.polymer_model, 
        polymer_phase_density, 
        temperature, 
        polymer_phase_mass_fractions)
    return a
end

function solve_polymer_density(
    model::DGRPTModel, temperature::Number, polymer_phase_mass_fractions::AbstractVector; 
    initial_density=nothing, taylor_series_order=default_dgrpt_taylor_expansion_order, method=:roots)

    max_polymer_density = density_upper_bound(model.polymer_model, polymer_phase_mass_fractions) * polymer_phase_mass_fractions[1]

    if isnothing(initial_density)
        initial_density = model.polymer_dry_density * polymer_phase_mass_fractions[1]

        # sometimes the dry density is more dense than the theoretical maximum density of the polymer at that composition (e.g.,, in most liquids and swelling vapors)
        if initial_density >= max_polymer_density  
            initial_density = (max_polymer_density * 0.999)
        end
    end

    if method==:roots
        # roots solution 1
        roots_polymer_density_target = make_roots_polymer_density_target(model, temperature, polymer_phase_mass_fractions; taylor_series_order)
        
        # if typeof(max_polymer_density) <: ForwardDiff.Dual
        #     @show max_polymer_density.value
        #     @show isnan(max_polymer_density)
        # else
        #     @show max_polymer_density
        # end
        # if eltype(polymer_phase_mass_fractions) <: ForwardDiff.Dual
        #     @show [polymer_phase_mass_fraction.value for polymer_phase_mass_fraction in polymer_phase_mass_fractions]
        # else
        #     @show max_polymer_density
        # end
        
        solved_density = find_zeros(roots_polymer_density_target, eps()*10, max_polymer_density - eps(); no_pts=12, naive=false)
        if length(solved_density) > 1
            # @warn string(length(solved_density)) * " densities found! Choosing the maximum of the two." # todo deal with this later
            # if eltype(solved_density) <: ForwardDiff.Dual
            #     @show [s.value for s in solved_density], [root.value for root in roots_polymer_density_target.(solved_density)] 
            # else
            #     @show solved_density
            # end

            # if eltype(polymer_phase_mass_fractions) <: ForwardDiff.Dual
            #     @show [s.value for s in polymer_phase_mass_fractions]
            # else
            #     @show polymer_phase_mass_fractions
            # end
            return maximum(solved_density)
        
        elseif length(solved_density) == 0
            # @warn "No densities found! Defaulting to the Optim method"
            # if eltype(polymer_phase_mass_fractions) <: ForwardDiff.Dual
            #     @show [s.value for s in polymer_phase_mass_fractions]
            # else
            #     @show polymer_phase_mass_fractions
            # end
            return solve_polymer_density(model, temperature, polymer_phase_mass_fractions; taylor_series_order, method=:optim)

        else
            return solved_density[1]
        end
    elseif method==:optim
        #optim solution 1
        optim_polymer_density_target = make_optim_polymer_density_target(model, temperature, polymer_phase_mass_fractions; taylor_series_order)
        # @show max_polymer_density
        res = Optim.optimize(
            optim_polymer_density_target, [eps()], [max_polymer_density -  eps()],
            [initial_density], 
            Fminbox(LBFGS()), 
            Optim.Options(
                allow_f_increases = false,
                g_tol = 1e-8
            ); autodiff=:forward)
        solved_density = Optim.minimizer(res)[1]
        return solved_density
    end
end

function expected_polymer_chemical_potential(model::DGRPTModel, temperature, polymer_density, polymer_phase_mass_fractions; expansion_order=default_dgrpt_taylor_expansion_order)
    # component_chemical_potential_expansions = Vector(zeros(length(model.penetrants)))  # can't use staticvectors with forwarddiff :(
    num_penetrants = length(polymer_phase_mass_fractions) - 1

    component_chemical_potential_expansions = 0

    for pen_idx in eachindex(polymer_phase_mass_fractions[2:end])
        function polymer_chemical_potential_taylor_expansion(penetrant_mass_fraction)
            mass_fracs = vcat([1 - penetrant_mass_fraction], zeros(num_penetrants))
            mass_fracs[pen_idx + 1] = penetrant_mass_fraction 
            chem_pot = calculate_polymer_phase_chemical_potentials(
                model, 
                temperature,
                polymer_density, 
                mass_fracs)[1]
            return chem_pot
        end

        penetrant_taylor_expansion_function = expand_taylor_higher_terms(
            polymer_chemical_potential_taylor_expansion,
            1e-6, expansion_order
        )
        evaluated_penetrant_taylor_expansion = penetrant_taylor_expansion_function(polymer_phase_mass_fractions[pen_idx + 1])
        component_chemical_potential_expansions += evaluated_penetrant_taylor_expansion
    end
    component_chemical_potential_expansions
    dry_μ = dry_polymer_chemical_potential(model, temperature, length(polymer_phase_mass_fractions))

    return dry_μ + sum(component_chemical_potential_expansions)
end

function dry_polymer_chemical_potential(model::DGRPTModel, temperature, num_components)
    pseudo_mass_fracs = zeros(Float64, num_components)
    pseudo_mass_fracs[1] = 1
    μ = ρTω_chemical_potential(
        model.polymer_model, 
        model.polymer_dry_density, 
        temperature, 
        pseudo_mass_fracs)
    return μ[1]
end

# Target functions for optim and roots

function make_penetrant_mass_fraction_target(
    model::DGRPTModel, 
    temperature::Number, 
    pressure::Number, 
    bulk_penetrant_mole_fractions::AbstractVector{<:Number}; 
    taylor_series_order=default_dgrpt_taylor_expansion_order)
    
    target_penetrant_activities = calculate_bulk_phase_activities(model, temperature, pressure, bulk_penetrant_mole_fractions)
    activity_scaling_factors = 1 ./ target_penetrant_activities

    target = function error_function(penetrant_mass_fractions)
        polymer_mass_fraction = 1 - sum(penetrant_mass_fractions)
        polymer_phase_mass_fractions = vcat(polymer_mass_fraction, penetrant_mass_fractions)
        polymer_density = solve_polymer_density(model, temperature, polymer_phase_mass_fractions; taylor_series_order)
        polymer_phase_activities = calculate_polymer_phase_activities(model, temperature, polymer_density, polymer_phase_mass_fractions)
        residual_squared = rss(target_penetrant_activities .* activity_scaling_factors, polymer_phase_activities[2:end] .* activity_scaling_factors)    
        
        # if eltype(penetrant_mass_fractions) <: ForwardDiff.Dual
        #     @show [frac.value for frac in penetrant_mass_fractions]
        # else
        #     @show penetrant_mass_fractions
        # end
        # if typeof(residual_squared) <: ForwardDiff.Dual
        #     @show residual_squared.value
        # else
        #     @show residual_squared
        # end
        
        return residual_squared 
    end

    return target
end

function make_roots_polymer_density_target(model::DGRPTModel, temperature::Number, polymer_phase_mass_fractions::AbstractVector; 
    taylor_series_order=default_dgrpt_taylor_expansion_order)
    
    scaling_factor = 1 / expected_polymer_chemical_potential(
        model, temperature, model.polymer_dry_density, polymer_phase_mass_fractions; expansion_order=taylor_series_order)
    
    function roots_polymer_density_target(density)  # target for roots
        target_polymer_chemical_potential = expected_polymer_chemical_potential(
            model, temperature, density, polymer_phase_mass_fractions; expansion_order=taylor_series_order)

        polymer_chemical_potential = PolymerMembranes.calculate_polymer_phase_chemical_potentials(
            model, temperature, density, polymer_phase_mass_fractions)[1]
        return (target_polymer_chemical_potential - polymer_chemical_potential) * scaling_factor
    end
    return roots_polymer_density_target
end

function make_optim_polymer_density_target(model::DGRPTModel, temperature::Number, polymer_phase_mass_fractions::AbstractVector; 
    taylor_series_order=default_dgrpt_taylor_expansion_order)
    
    scaling_factor = 1 / expected_polymer_chemical_potential(
        model, temperature, model.polymer_dry_density, polymer_phase_mass_fractions; expansion_order=taylor_series_order)

    function optim_polymer_density_target(density)  # target for optim
        target_polymer_chemical_potential = expected_polymer_chemical_potential(
            model, temperature, density[1], polymer_phase_mass_fractions; expansion_order=taylor_series_order)

        polymer_chemical_potential = PolymerMembranes.calculate_polymer_phase_chemical_potentials(
            model, temperature, density[1], polymer_phase_mass_fractions)[1]
        return ((target_polymer_chemical_potential - polymer_chemical_potential) * scaling_factor)^2
    end
    return optim_polymer_density_target
end
