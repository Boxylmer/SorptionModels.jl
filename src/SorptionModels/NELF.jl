const DEFAULT_NELF_INFINITE_DILUTION_PRESSURE = 1e-8 # ???
const DEFAULT_NELF_POLYMER_MOLECULAR_WEIGHT = 100000

struct NELF end
# struct NELF{M} end
# NELF() = NELF(:unspecified) # TODO

"""
    NELFModel(bulk_model, polymer_model, polymer_dry_density)

Create a NELF sorption model, using two equations of state, one representing the bulk phase and one representing the polymer phase. 

Notes:
- The `polymer_dry_density` should reflect the density of the pure polymer at STP. 
- Requires that these models contain the polymer is the first index when referencing compositions and interaction parameters. The remaining indexes must match `penetrants` in the predictions functions.

`F. Doghieri, G.C. Sarti, Nonequilibrium Lattice Fluids: A Predictive Model for the Solubility in Glassy Polymers, Macromolecules. 29 (1996) 7885–7896. https://doi.org/10.1021/ma951366c.`
"""
struct NELFModel{BMT, POLYMT, PDT} <: SorptionModel
    bulk_model::BMT                 # EOS
    polymer_model::POLYMT           # EOS
    polymer_dry_density::PDT        # number
end

function NELFModel(polymer_phase_model, polymer_dry_density)
    bulk_phase_model = Clapeyron.split_model(polymer_phase_model, [2:length(polymer_phase_model)])[1]
    return NELFModel(bulk_phase_model, polymer_phase_model, polymer_dry_density)
end

function predict_concentration(model::NELFModel, temperature::Number, pressure::Number, bulk_phase_mole_fractions=[1]; ksw=nothing, units=:cc)
    minimum_val = 100 * eps()
    error_target = _make_nelf_model_mass_fraction_target(model, temperature, pressure, bulk_phase_mole_fractions; ksw, minimum_val)

    penetrant_mass_fraction_initial_guesses = zeros(length(bulk_phase_mole_fractions)) .+ 0.0001
    lower = zeros(length(penetrant_mass_fraction_initial_guesses))
    upper = ones(length(penetrant_mass_fraction_initial_guesses)) .- eps()

    res = Optim.optimize(
        error_target, 
        lower, upper,
        penetrant_mass_fraction_initial_guesses, 
        Fminbox(NelderMead()),
        Optim.Options(
        allow_f_increases = false,
        ); autodiff=:forward
    )

    penetrant_mass_fractions = Optim.minimizer(res)
    polymer_phase_mass_fractions = vcat(1 - sum(penetrant_mass_fractions), penetrant_mass_fractions)

    if units==:frac
        return polymer_phase_mass_fractions[2:end]
    elseif units==:g
        concs_g_g = polymer_phase_mass_fractions_to_gpen_per_gpol(polymer_phase_mass_fractions)
        return concs_g_g
    elseif units==:cc
        concs_cc_cc = polymer_phase_mass_fractions_to_ccpen_per_ccpol(polymer_phase_mass_fractions, model.polymer_dry_density, Clapeyron.mw(model.bulk_model))
        return concs_cc_cc
    else
        throw(ArgumentError("Units not supported: " * string(units)))
    end
end

function _make_nelf_model_mass_fraction_target(model::NELFModel, temperature::Number, pressure::Number, bulk_phase_mole_fractions; ksw=nothing, minimum_val=100*eps())
    if isnothing(ksw)
        ksw = zeros(length(bulk_phase_mole_fractions))
    end

    pressure = max(minimum_val * one(pressure), pressure)    
    bulk_phase_mole_fractions = ((i) -> (i < minimum_val ? minimum_val : i)).(bulk_phase_mole_fractions)  # Courtesy of Clementine (Julia Discord)

    target_potentials = Clapeyron.chemical_potential(model.bulk_model, pressure * MembraneBase.PA_PER_MPA, temperature, bulk_phase_mole_fractions)
   
    normalizer = (rss(zeros(length(bulk_phase_mole_fractions)), target_potentials))

    function error_target(penetrant_mass_fractions)
        polymer_mass_fraction = 1 - sum(penetrant_mass_fractions)

        polymer_phase_mass_fractions = vcat(polymer_mass_fraction, penetrant_mass_fractions)
        polymer_phase_density_after_swelling = calculate_polymer_phase_density(model, pressure, bulk_phase_mole_fractions, polymer_phase_mass_fractions, ksw)
        polymer_phase_mole_fractions = mass_fractions_to_mole_fractions(polymer_phase_mass_fractions, Clapeyron.mw(model.polymer_model))

        # you have to calculate this every time. There is no other option when doing multicomponent as the phase maximum density changes at any composition.
        
        # TODO don't bother calculating mole frac and just calculate based on mass frac
        # polymer_phase_density_upper_bound = ub_density(model.polymer_model, polymer_phase_mole_fractions, :molar) #g/cm3

        # if polymer_phase_density_after_swelling > polymer_phase_density_upper_bound
        #     # println(".")
        #     # return NaN 
        #     return log1p(normalizer * polymer_phase_density_after_swelling) 
        # end

        # upper_bound_penalty = log(1/(polymer_phase_density_upper_bound - polymer_phase_density_after_swelling))

        polymer_phase_residual_potential = ρTw_chemical_potential( 
            model.polymer_model, 
            polymer_phase_density_after_swelling, 
            temperature, 
            polymer_phase_mole_fractions,
            :molar)
        residual_squared = log1p((rss(target_potentials, polymer_phase_residual_potential[2:end])))
        # println(residual_squared)
        return residual_squared #* upper_bound_penalty
    end
    return error_target
end


function calculate_swelled_polymer_density(model::NELFModel, penetrant_partial_pressures::AbstractVector{<:Number}, ksw_values::AbstractVector{<:Number})
    return model.polymer_dry_density * (1 - sum(ksw_values .* penetrant_partial_pressures))
end

function calculate_swelled_polymer_density(model::NELFModel, pressure::Number, penetrant_mole_fractions::AbstractVector{<:Number}, ksw_values::AbstractVector{<:Number})
    partial_pressures = penetrant_mole_fractions .* pressure
    return calculate_swelled_polymer_density(model, partial_pressures, ksw_values)
end

function calculate_polymer_phase_density(swelled_density::Number, polymer_phase_mass_fractions::AbstractVector{<:Number})
    return swelled_density / polymer_phase_mass_fractions[1]
end

function calculate_polymer_phase_density(model::NELFModel, pressure::Number, bulk_phase_mole_fractions::AbstractVector{<:Number}, polymer_phase_mass_fractions::AbstractVector{<:Number}, ksw_values::AbstractVector{<:Number})
    polymer_density_after_swelling = calculate_swelled_polymer_density(model, pressure, bulk_phase_mole_fractions, ksw_values)
    return calculate_polymer_phase_density(polymer_density_after_swelling, polymer_phase_mass_fractions)
end

"""
    infinite_dilution_solubility(model::NELFModel, temperature::Number)
Get infinite dilution solubility in **((CC/CC) / MPa)**
"""
# todo how will we predict the infinite dilution solubility or reject multicomponent models?
function infinite_dilution_solubility(
    model::NELFModel, temperature::Number; 
    infinite_dilution_pressure = DEFAULT_NELF_INFINITE_DILUTION_PRESSURE)  # naieve
    return predict_concentration(model, temperature, infinite_dilution_pressure, [1]; ksw=[0])[1] / infinite_dilution_pressure
end


# functions for fit data to NELF parameters
"""
    fit_model(::NELF, model_choice, isotherms, bulk_phase_characteristic_params, [polymer_molecular_weight]; [custom_densities], [initial_search_resolution], [uncertainty_method])
Find the EOS parameters of a polymer from a vector of `IsothermData`s using the NELF model. 
# Arguments
- `model_choice`: MembraneEOS model to use
- `isotherms`: `Vector` of `IsothermData` structs using the polymer in question. Temperature, pressure, density, and concentration must be provided.
- `bulk_phase_params`: Vector of pure characteristic parameter vectors following the same order as the isotherms. 
  - Each parameter vector's units should be formatted as (P★, T★, ρ★, mw) in units of (mpa, K, g/cm3, g/mol).
  - E.g., for Sanchez Lacombe and two input isotherms, `bulk_phase_characteristic_params = [[p★_1, t★_1, ρ★_1, mw_1], [p★_2, t★_2, ρ★_2, mw_2]]`
- `polymer_molecular_weight`: A known molecular weight for the polymer (otherwise a default, arbitrarily large value is assumed)
- `custom_densities`: An array (matching the dimensions of the input isotherms) of densities to use instead of the ones in the isotherm data provided. 
- `uncertainty_method`: Calculate the uncertainty of the parameters from the fitting. For NELF, the Hessian method is currently the only one implemented.
- `kij_fitting_method`: :infinite_dilution, :all, :none
"""
function fit_model(::NELF,# TODO Docs update
    isotherms::AbstractVector{<:IsothermData}, 
    bulk_phase_params::AbstractVector;
    polymer_molecular_weight=DEFAULT_NELF_POLYMER_MOLECULAR_WEIGHT, 
    verbose=false, 
    initial_search_resolution=20,
    custom_densities::Union{Missing, AbstractArray}=missing, 
    uncertainty_method=nothing, 
    adjust_kij=false)
    
    if verbose
        println("Starting parameter generation for NELF fit with the Sanchez Lacombe EoS")
    end
    # this function uses SL, which needs 4 params per component, one of which is already specified (MW)

    error_function = _make_nelf_model_parameter_target(
        isotherms, 
        bulk_phase_params, 
        DEFAULT_NELF_INFINITE_DILUTION_PRESSURE, 
        polymer_molecular_weight; 
        adjust_kij)

    preliminary_error_function = _make_nelf_model_parameter_target(
        isotherms, 
        bulk_phase_params, 
        DEFAULT_NELF_INFINITE_DILUTION_PRESSURE, 
        polymer_molecular_weight;
        adjust_kij=false)

    if ismissing(custom_densities)
        densities = polymer_density.(isotherms)
    else
        densities = custom_densities
    end
    density_lower_bound = maximum(densities)
    density_upper_bound = 3
    naive_lower = [50., 50., density_lower_bound]
    naive_upper = [2000., 2000., density_upper_bound]
    if verbose
        println("Identified naive lower bound of $naive_lower and upper bound of $naive_upper for P, T and rho respectively.")
    end
    min_results, min_args, lower, upper = scan_for_starting_point_and_bounds_3_dims(preliminary_error_function, naive_lower, naive_upper, initial_search_resolution; verbose)
    if verbose
        println("Found initial starting point at point $min_args, with an RSS of $min_results, bounded by a lower bound of $lower and upper bound of $upper. Starting optimization...")
    end
    # @show min_args, lower, upper
    res = Optim.optimize(
        error_function, 
        # lower, upper, 
        # [500, 500, density_lower_bound * 1.2], 
        min_args,
        NelderMead(), 
        # Fminbox(LBFGS()), 
        Optim.Options(; allow_f_increases = false))

    _minimizer = Optim.minimizer(res)

    if uncertainty_method == :jackKnife
        throw(ErrorException("Jackknifing isn't working for this set of models just yet"))
        errors = jackknife_uncertainty(fitting_uncertainty_wrapper, dataset(_step_data))

    elseif uncertainty_method == :Hessian
        errors = rss_minimizer_standard_errors(error_function, _minimizer, length(isotherms))
        minimizer = [_minimizer[i] .± errors[i] for i in eachindex(_minimizer, errors)]
    elseif isnothing(uncertainty_method)
        minimizer = _minimizer
    else
        throw(ArgumentError("Uncertainty method should be a symbol, e.g., :Hessian, or :JackKnife"))
    end
    return [minimizer..., polymer_molecular_weight]
    # work in progress
end

# TODO possibly move to model methods
"Construct a two component SL EoSModel with params formatted as (P★, T★, ρ★, mw) in units of (mpa, K, g/cm3, g/mol)."
function construct_binary_sl_eosmodel(polymer_params, bulk_params, kijval::Number)
    v★(P★, T★,) = 8.31446261815324 * T★ / P★ / 1000000 # J / (mol*K) * K / mpa -> pa * m3 / (mol * mpa) ->  need to divide by 1000000 to get m3/mol
    ϵ★(T★) = 8.31446261815324 * T★ # J / (mol * K) * K -> J/mol 
    r(P★, T★, ρ★, mw) = mw * (P★ * 1000000) / (8.31446261815324 * T★ * (ρ★ / 1e-6)) # g/mol * mpa * 1000000 pa/mpa / ((j/mol*K) * K * g/(cm3 / 1e-6 m3/cm3)) -> unitless

    P★, T★, ρ★, mw = [[pol, bulk] for (pol, bulk) in zip(polymer_params, bulk_params)]
    # @show P★, T★, ρ★, mw = collect(zip(polymer_params, bulk_params))
    kij_mat = [0 kijval; kijval 0]
    model = Clapeyron.SL(
        ["Polymer", "Component"], 
        userlocations = Dict(
            :vol => v★.(P★, T★), 
            :segment => r.(P★, T★, ρ★, mw),
            :epsilon => ϵ★.(T★), 
            :Mw => mw,
            # :k => [0 kijval; kijval 0]
            
        ),
        mixing_userlocations = (;k0 = kij_mat, k1 = [0 0; 0 0], l = [0 0; 0 0])
    )
    return model
end

# todo this allocates like crazy when it really shouldn't. What's going on?
"""
    scan_for_starting_point_and_bounds_3_dims(target_function::Function, naive_lower::Vector{Float64}, naive_lower::Vector{Float64}, resolutions=missing; return_grid=false, verbose=true)

    Search for a vector of parameters, bounded by `naive_lower` and `naive_lower`, that is closest to minimizing a `target_function` by trying every possible value in a grid of `resolutions`. 
- `resolutions` can be a vector of the same length as the bounds and arguments to the `target_function`. If a single integer is passed, it will assume you want the same resolution on all input dimensions.
- This search algorithm assumes that the target_function contains one obvious local minima, but is robust to `NaN`, `missing`, and `nothing` output from `target_function`.

Returns `min_results`, `min_args`, `lower_bounds`, `upper_bounds`, or rather, the output of the function that was searched, the argument vector to get said output, and the lower and upper bounds in which a better solution may exist.
Returns only the error grid that was evaluated and skips finding any other results if `return_grid` is true. (used for diagnostic purposes)

    Example
```
function expensive_func(x)
    return (x[1]-15)^2 + (x[2]-10)^2 + (x[3]-12)^2
end
@show scan_for_starting_point_and_bounds_3_dims(expensive_func, (1.0, 1.0, 1.), (20., 20., 20.), (10, 11, 12); verbose=true)
```
# """
function scan_for_starting_point_and_bounds_3_dims(target_function::Function, naive_lower, naive_upper, resolutions=20; verbose=false)
    @assert length(naive_lower) == length(naive_upper)
    dim_resolutions = MVector{length(naive_lower), Int64}(undef)
    if typeof(resolutions) <: Number
        dim_resolutions .= resolutions, resolutions, resolutions
    else
        dim_resolutions .= resolutions
    end
    n_dims = length(dim_resolutions)
    
    get_arg_value(dimension::Int64, index::Int64)::Float64 = naive_lower[dimension] + (naive_upper[dimension]-naive_lower[dimension])*(index/dim_resolutions[dimension])
    function get_args(indices)
        mvec = MVector{length(indices), Float64}(undef)
        get_args!(mvec, indices)
        return mvec
    end
    function get_args!(mvec, indices)
        for idx in eachindex(indices, mvec)
            mvec[idx] = get_arg_value(idx, indices[idx])
        end
    end
    
    # write in what should get allocated
    min_results = Inf                               # no allocation somehow
    min_indices = MVector{n_dims, Int64}(undef)     # 1 allocation
    temp_args = MVector{n_dims, Float64}(undef)     # 1 allocation
    temp_indices = MVector{n_dims, Int64}(undef)    # 1 allocation
    completed_iters = 0
    
    for p_index in 1:dim_resolutions[1]
        for t_index in 1:dim_resolutions[2]
            for ρ_index in 1:dim_resolutions[3]
                temp_indices .= p_index, t_index, ρ_index
                get_args!(temp_args, temp_indices)
                res = target_function(temp_args)

                if verbose; completed_iters += 1; end
                if isnothing(res) || ismissing(res) || isnan(res); continue; end
                if res < min_results
                    min_results = res
                    min_indices .= temp_indices
                end
                if verbose && (mod(completed_iters, 500) == 0)
                    print("Initial naive search is " * string(round(completed_iters / (dim_resolutions[1] * dim_resolutions[2] * dim_resolutions[3]) * 100))*"% complete.\r")
                end
            end
        end
    end
    min_args = get_args(min_indices)                    # 1 allocation
    min_indices .+= 1
    upper_bounds = get_args(min_indices)                # 1 allocation
    min_indices .-= 2
    lower_bounds = get_args(min_indices)                # 1 allocation
    return min_results, min_args, lower_bounds, upper_bounds
end

function _make_nelf_model_parameter_target(
    isotherms::AbstractVector{<:IsothermData}, 
    bulk_phase_param_vecs::AbstractVector, 
    infinite_dilution_pressure=DEFAULT_NELF_INFINITE_DILUTION_PRESSURE, 
    polymer_molecular_weight=100000;
    adjust_kij=false) 
    dualmode_models = [fit_model(DualMode(), isotherm) for isotherm in isotherms] 
    # TODO don't restrict to dual mode type isotherms only, we probably want to just pick the first data point maybe? Maybe need a flag to define strategy
    
    densities = polymer_density.(isotherms) # get each isotherm's density in case the user accounted for polymers from different batches
    temperatures = temperature.(isotherms)

    
    if !adjust_kij # TODO the strategy below would make this faster (just don't fit kij)
        function error_function_ideal_kij(char_param_vec::AbstractVector{<:Number})
            given_val = zeros(length(isotherms))
            pred_val = zeros(length(isotherms)) 
            for i in eachindex(isotherms, bulk_phase_param_vecs)
                kij = 0
                polymer_phase_model = construct_binary_sl_eosmodel([char_param_vec..., polymer_molecular_weight], bulk_phase_param_vecs[i], kij) 
    
                nelf_model = NELFModel(polymer_phase_model, densities[i]) 
                pred_val[i] = predict_concentration(nelf_model, temperatures[i], infinite_dilution_pressure)[1]
                given_val[i] = predict_concentration(dualmode_models[i]::DualModeModel, infinite_dilution_pressure) 
                if isnan(pred_val[i])
                    println("Nan found")
                end
               
                # working version
                # pred_val[i] = predict_concentration(nelf_model, temperatures[i], infinite_dilution_pressure, [1]; ksw=[0], nan_on_failure)[1] #/ infinite_dilution_pressure
                # given_val[i] = predict_concentration(dualmode_models[i]::DualModeModel, infinite_dilution_pressure)
                
               
            end
            resid = sum(((given_val .- pred_val) ./ given_val).^2)
            err = log1p(resid)
            return err
        end
        return error_function_ideal_kij
    else
        index_map = map_equivalent_component_indices(bulk_phase_param_vecs)
        # current_kij_values = Vector{Float64}(undef, length(index_map))
        function error_function_fit_kij(char_param_vec::AbstractVector{<:Number})
            given_val = zeros(length(isotherms))
            pred_val = zeros(length(isotherms)) 
            for (component_parameters, param_indices) in index_map
                relevant_isotherms = @view isotherms[param_indices]
                relevant_interpolation_models = @view dualmode_models[param_indices]
                eosmodel = construct_binary_sl_eosmodel([char_param_vec..., polymer_molecular_weight], component_parameters, 0.0)
                nelf_model = NELFModel(eosmodel, densities[1]) # TODO this is just taking densities at point one, need to fix in future
                fit_kij!(nelf_model, relevant_isotherms, relevant_interpolation_models; infinite_dilution_pressure)

                for parameter_idx in param_indices
                    isotherm = isotherms[parameter_idx]
                    p = MembraneBase.pressures(isotherm; component=1)[findfirst(x -> x != 0, MembraneBase.pressures(isotherm; component=1))]
                    pred_val[parameter_idx] = predict_concentration(nelf_model, temperatures[parameter_idx], p)[1] 
                    given_val[parameter_idx] = predict_concentration(dualmode_models[parameter_idx]::DualModeModel, p) 
                end
            end
            resid = sum(((given_val .- pred_val) ./ given_val).^2)
            err = log1p(resid)
            return err
        end
        return error_function_fit_kij
    end
end

# TODO Docs
""" 
    fit_kij(NELF(), isotherms, bulk_parameters, polymer_parameters; [interpolation_model]=DualMode(), [kij_fit_p_mpa]=1e-4)
Find the best kij parameter for one component in the NELF model.  
# Arguments
- `isotherms`: `Vector` of `IsothermData` structs using the polymer in question with **only one kind of component**. Temperature, pressure, density, and concentration must be provided.
- `bulk_parameters`: vector of the gas or vapor's characteristic parameters, in the format of `[p★_mpa, t★_k, ρ★_g_cm3, mw_g_mol]`.
- `polymer_parameters`: vector of the polymer's characteristic parameters, in the format of `[p★_mpa, t★_k, ρ★_g_cm3, mw_g_mol]`.
"""
function fit_kij(::NELF, isotherms::AbstractArray{<:IsothermData}, 
        polymer_parameters::AbstractVector{<:Number},
        bulk_parameters::AbstractVector{<:Number}; 
        kwargs...
    )
    
    polymer_phase_model = construct_binary_sl_eosmodel(polymer_parameters, bulk_parameters, 0.0)
    densities = polymer_density.(isotherms) # get each isotherm's density in case the user accounted for polymers from different batches
    # right now just warn the user if they're different
    if !all((densities[1] == densities[i] for i in eachindex(1:length(densities))))
        throw(ArgumentError("Isotherms had variable densities, this isn't currently handled!")) 
    end
    
    nelf_model = NELFModel(polymer_phase_model, densities[1])
    fit_kij!(nelf_model, isotherms; kwargs...)
end

function fit_kij!(nelfmodel::NELFModel,
        isotherms::AbstractArray{<:IsothermData}; 
        interpolation_model=DualMode(), 
        infinite_dilution_pressure=DEFAULT_NELF_INFINITE_DILUTION_PRESSURE
    )

    interpolation_models = [fit_model(interpolation_model, isotherm) for isotherm in isotherms]
   return fit_kij!(nelfmodel, isotherms, interpolation_models; infinite_dilution_pressure)
end

function fit_kij!(nelfmodel::NELFModel,
        isotherms::AbstractArray{<:IsothermData},
        interpolation_models::AbstractArray{<:SorptionModel}; 
        infinite_dilution_pressure=DEFAULT_NELF_INFINITE_DILUTION_PRESSURE
    )
    temperatures = temperature.(isotherms)

    function kij_error_function(kij)
        kij_mat = [0 kij[1]; kij[1] 0]
        Clapeyron.set_k!(nelfmodel.polymer_model, kij_mat)
        pred_concs = [predict_concentration(nelfmodel, temperatures[i], infinite_dilution_pressure)[1] for i in eachindex(interpolation_models)] ./ infinite_dilution_pressure
        exp_concs = predict_concentration.(interpolation_models, infinite_dilution_pressure) ./ infinite_dilution_pressure
        err = rss(exp_concs, pred_concs)
        return err
    end
    return Optim.minimizer(Optim.optimize(kij_error_function, [0.], NelderMead()))[1]
end

function fit_ksw(::NELF, isotherm::IsothermData, polymer_model, bulk_model)
    density = polymer_density(isotherm)
    temp = temperature(isotherm)
    pressures_mpa = partial_pressures(isotherm; component=1)
    concentrations_cc = concentration(isotherm; component=1)
    nelf_model = NELFModel(bulk_model, polymer_model, density)
    function error_function(ksw)
        pred_concs = [predict_concentration(nelf_model, temp, p, [1]; ksw=ksw)[1] for p in pressures_mpa]
        @show ksw
        return rss(concentrations_cc, pred_concs)
    end
    # lower = [0.]
    # upper = [Inf]
    init = [0.0]
    @show res = Optim.optimize(error_function, init, BFGS())
    @show minimizer = Optim.minimizer(res)
    return min(0.0, minimizer)
end


"Given a vector of eos parameters, return a set of vectors that contain the indices of identical parameters."
function map_equivalent_component_indices(params::AbstractVector{<:AbstractVector{<:Number}})
    index_map = Dict{Tuple, Vector{Int64}}()
    for (i, param_vector) in enumerate(params)
        param_tuple = tuple(param_vector...)
        if !haskey(index_map, param_tuple)
            index_map[param_tuple] = Vector{Int64}()
        end 
        push!(index_map[param_tuple], i)
    end
    index_map
end

