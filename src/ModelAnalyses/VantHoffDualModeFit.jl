struct VantHoffDualModeAnalysis  # todo parameterization, this is really just a container that's called once (when optimized)
    ch0_final
    b0_final
    kd0_final
    mch_final
    ΔHb_final
    ΔHkd_final
    ch0_initial
    b0_initial
    kd0_initial
    mch_initial
    ΔHb_initial
    ΔHkd_initial
    given_isotherms
    final_predicted_isotherms
    initial_predicted_isotherms
    final_rss
    initial_rss
    final_models
    initial_models
    standalone_models
    covariance_matrix
    use_fugacity::Bool
end

module VHDMFHelper 
    using PolymerMembranes 
    using ..PolymerMembranes: R_J_MOL_K, rss, strip_measurement_to_value, fit_linear_data, IsothermData, TPCDataset, get_isotherms, jackknife_uncertainty
    using ..PolymerMembranes: rss_covariance_matrix, rss_minimizer_standard_errors, approximate_hessian, inverse_hessian
    using StaticArrays
    using Measurements
    using Optim
    
    # param_vec = [ch0, b0, kd0, mch, ΔHb, ΔHkd]     
    function get_model_vectors(param_vec::AbstractVector, temps::AbstractVector{<:Number})
        ch_values = SVector{length(temps), eltype(param_vec)}([param_vec[1] + param_vec[4] * temp for temp in temps]...)
        b_values = SVector{length(temps), eltype(param_vec)}([param_vec[2] * exp(-param_vec[5] / (R_J_MOL_K * temp)) for temp in temps]...)
        kd_values = SVector{length(temps), eltype(param_vec)}([param_vec[3] * exp(-param_vec[6] / (R_J_MOL_K * temp)) for temp in temps]...)
        model_vectors = SVector{length(temps), Vector{eltype(param_vec)}}([[ch_values[idx], b_values[idx], kd_values[idx]] for idx in eachindex(ch_values, b_values, kd_values)]...)
        return model_vectors 
    end

    function initialize_linearized_parameter_vector(isotherms::AbstractVector{<:IsothermData}; use_fugacity=false)
        # fit the initial models and use them as the initial parameters for the general optimization
        standalone_models = [fit_model(DualMode(), isotherm; use_fugacity=use_fugacity) for isotherm in isotherms]
        chs = [initial_model.ch for initial_model in standalone_models] 
        bs = [initial_model.b for initial_model in standalone_models] 
        kds = [initial_model.kd for initial_model in standalone_models] 
        
        # make some pre-conditioned guesses for the optimization
        temps = strip_measurement_to_value([temperature(isotherm) for isotherm in isotherms])
        inverse_temps = 1 ./ temps
        log_bs = log.(bs)
        log_kds = log.(kds)
        negative_ΔHb_per_R, log_b0 = strip_measurement_to_value.(fit_linear_data(inverse_temps, log_bs))
        negative_ΔHkd_per_R, log_kd0 = strip_measurement_to_value.(fit_linear_data(inverse_temps, log_kds))
        ΔHb = min(-10000, -negative_ΔHb_per_R * R_J_MOL_K) / 1000 
        ΔHkd = min(-10000, -negative_ΔHkd_per_R * R_J_MOL_K) / 1000
        log_b0 = max(-10, min(-4., log_b0))  # for SAMIN()
        log_kd0 = max(-10, min(-4., log_kd0)) # for SAMIN()
        ch0 = max(chs...)
        mch = -0.01
        linearized_initial_param_vec = [ch0, log_b0, log_kd0, mch, ΔHb, ΔHkd]
        # initial_param_vec = unlinearize_input_vector(linearized_initial_param_vec)  # save for later
        return linearized_initial_param_vec
    end

    function linearize_input_vector(vec::AbstractVector)
        return @SVector [
            vec[1], log(vec[2]), log(vec[3]), vec[4], vec[5] / 1000, vec[6] / 1000
        ]
    end
    
    function unlinearize_input_vector(vec::AbstractVector)
        return @SVector [
            vec[1], exp(vec[2]), exp(vec[3]), vec[4], vec[5] * 1000, vec[6] * 1000
        ]
    end

    # rss given a model vector
    function total_rss(model_vectors::AbstractVector{<:AbstractVector}, isotherms, use_fugacity::Bool)
        models = SVector{length(model_vectors), DualModeModel}([DualModeModel(model_vector...; use_fugacity) for model_vector in model_vectors])
        err = sum([rss(models[idx], isotherms[idx]; use_fugacity=use_fugacity) for idx in eachindex(models, isotherms)])
        return (err)
    end

    # rss given a param vector guess
    function total_rss(param_vec::AbstractVector{<:Number}, isotherms, use_fugacity::Bool)
        temps = SVector{length(isotherms), Float64}(strip_measurement_to_value([temperature(isotherm) for isotherm in isotherms]))
        model_vectors = get_model_vectors(param_vec, temps)
        total_error = total_rss(model_vectors, isotherms, use_fugacity)
        return total_error
    end

    function fitted_params_and_covariance(isotherms::AbstractVector{<:IsothermData}; use_fugacity=false, linearized_initial_param_vec=nothing) 
        function target_via_linearized_params(linearized_param_vec::AbstractVector)
            unlinearized_param_vec = unlinearize_input_vector(linearized_param_vec)
            return total_rss(unlinearized_param_vec, isotherms, use_fugacity)
        end

        function target_via_real_params(param_vec::AbstractVector)
            return total_rss(param_vec, isotherms, use_fugacity)
        end

        function gradtarget!(∇, x)
            ∇[:] = Optim.NLSolversBase.gradient!(Optim.NLSolversBase.OnceDifferentiable(target_via_linearized_params, x, autodiff = :finite), x)
            # ∇[:] = Optim.NLSolversBase.gradient!(Optim.NLSolversBase.OnceDifferentiable(forward_target, x, autodiff = :forward), x)
            return nothing
        end
        
        if isnothing(linearized_initial_param_vec)
            linearized_initial_param_vec = initialize_linearized_parameter_vector(isotherms; use_fugacity) 
        end

        # general requirementes like bounds and targets
        lower = [0., -16., -16., -100., -100., -100.]
        upper = [1000., 0., 0., 0., 0., 0.]
    
        # optim
        res = Optim.optimize(target_via_linearized_params, gradtarget!, lower, upper, linearized_initial_param_vec, Fminbox(BFGS()))  # doesn't find N2 minima well at all
 
        # res = Optim.optimize(target, lower, upper, linearized_initial_param_vec, SAMIN(; rt=0.9, ns=10, nt=10), Optim.Options(iterations = 10^7))
        linearized_final_param_vec = Optim.minimizer(res)
        final_param_vec = unlinearize_input_vector(linearized_final_param_vec)

        num_observations = sum([num_steps(isotherm) for isotherm in isotherms])
        param_standard_errors = rss_minimizer_standard_errors(target_via_real_params, final_param_vec, num_observations)
        final_params_with_errors = final_param_vec .± param_standard_errors
        covariance_matrix = rss_covariance_matrix(target_via_real_params, final_param_vec, num_observations)

        return final_params_with_errors, covariance_matrix
    end

end

"""
    temperature_dependent_dual_mode_fitting(isotherms::AbstractVector{<:IsothermData}; use_fugacity=false)

Fit a set of isotherms (using the same gas and polymer) at varying temperatures to the Dual Mode model. 
If `use_fugacity` is true, the solver will fit to fugacities within the isotherms rather than pressures.

Returns a `VantHoffDualModeAnalysis` struct containing contextual information on the fittings. 
Notable variables within it are:
- `ch0_final`: Converged Ch0 value.
- `b0_final`: Converged b0 value.
- `kd0_final`: Converged kd0 value.
- `mch_final`: Converged mch value.
- `ΔHb_final`: Converged ΔHb value.
- `ΔHkd_final`: Converged ΔHkd value.
- `final_rss`: Total sum of squared weighted residuals after convergence.
- `final_models`: Vector of Dual Mode models corresponding to the original order of the isotherms given.
- `covariance_matrix`: Matrix of covariances between individual parameters above (indexed in the original order they were listed)
"""
function temperature_dependent_dual_mode_fitting(isotherms::AbstractVector{<:IsothermData}; use_fugacity=false) 
    # generate the extra data for the analysis object
    linearized_initial_param_vec = VHDMFHelper.initialize_linearized_parameter_vector(isotherms; use_fugacity)
    initial_param_vec = VHDMFHelper.unlinearize_input_vector(linearized_initial_param_vec)
    final_param_vec, covariance_matrix = VHDMFHelper.fitted_params_and_covariance(isotherms; use_fugacity,)

    temps = strip_measurement_to_value([temperature(isotherm) for isotherm in isotherms])
    initial_model_vectors = VHDMFHelper.get_model_vectors(initial_param_vec, temps)
    initial_models = [DualModeModel(model_vec...; use_fugacity=use_fugacity) for model_vec in initial_model_vectors]
    final_model_vectors = VHDMFHelper.get_model_vectors(final_param_vec, temps)
    final_models = [DualModeModel(model_vec...; use_fugacity=use_fugacity) for model_vec in final_model_vectors]
    
    # save initial models and vectors to see where we started in the optimization
    if use_fugacity
        given_pressures_vectors = [fugacities(isotherm; component=1) for isotherm in isotherms]
    else
        given_pressures_vectors = [pressures(isotherm; component=1) for isotherm in isotherms]
    end

    # todo include standlone predictions
    standalone_models = [fit_model(DualMode(), isotherm; use_fugacity=use_fugacity) for isotherm in isotherms]

    initial_predictions = [
        predict_concentration(initial_models[idx], given_pressures_vectors[idx]) for idx in eachindex(initial_models, given_pressures_vectors)
    ]
    final_predictions = [
        predict_concentration(final_models[idx], given_pressures_vectors[idx]) for idx in eachindex(final_models, given_pressures_vectors)
    ]
    if use_fugacity
        initial_predicted_isotherms = [
            IsothermData(; temperature_k=temps[idx], concentrations_cc=initial_predictions[idx], fugacities_mpa=given_pressures_vectors[idx])
            for idx in eachindex(temps, initial_predictions)
        ]
        final_predicted_isotherms = [
            IsothermData(; temperature_k=temps[idx], concentrations_cc=final_predictions[idx], fugacities_mpa=given_pressures_vectors[idx])
            for idx in eachindex(temps, final_predictions)
        ]
    else
        initial_predicted_isotherms = [
            IsothermData(; temperature_k=temps[idx], concentrations_cc=initial_predictions[idx], partial_pressures_mpa=given_pressures_vectors[idx])
            for idx in eachindex(temps, initial_predictions)
        ]
        final_predicted_isotherms = [
            IsothermData(; temperature_k=temps[idx], concentrations_cc=final_predictions[idx], partial_pressures_mpa=given_pressures_vectors[idx])
            for idx in eachindex(temps, final_predictions)
        ]
    end
    initial_rss = VHDMFHelper.total_rss(initial_param_vec, isotherms, use_fugacity)
    final_rss = VHDMFHelper.total_rss(final_param_vec, isotherms, use_fugacity)

    analysis = VantHoffDualModeAnalysis(
        final_param_vec..., initial_param_vec..., isotherms, 
        final_predicted_isotherms, initial_predicted_isotherms, 
        final_rss, initial_rss, final_models, initial_models, standalone_models,
        covariance_matrix, use_fugacity)
    
    return analysis
end


# """
#     temperature_dependent_dual_mode_fitting(isotherms; use_fugacity=false)
# 
# Fits multiple isotherms of the same penetrant/polymer combination to models for predicting``C_{H}^{\\`}``, ``b``, and ``k_{d}`` values, which in turn predict the isotherms themselves.

# """
# struct VantHoffDualModeModelConstraints{HBT, HKDT, BT, KDT}
#     ΔHb::HBT
#     ΔHkd::HKDT
#     b0::BT
#     kd0::KDT
# end
# struct MultiGasVantHoffDualModeModelConstraints{MT, NBT, NKDT, HBT, HKDT}
#     M::MT
#     Nb::NBT
#     Nkd::NKDT
#     ΔHbs::HBT
#     ΔHkds::HKDT
# end

# module VHDMHelper
#     using ..PolymerMembranes
#     function get_b0(M, Nb, critical_temperature)
#         b0 = exp(Nb + M * critical_temperature)
#         return b0
#     end
#     function get_kd0(M, Nkd, critical_temperature)
#         kd0 = exp(Nkd + M * critical_temperature)
#         return kd0
#     end
#     function get_b(b0, ΔHb, temperature)
#         return b0 * exp(-ΔHb / (PolymerMembranes.R_J_MOL_K * temperature))
#     end
#     function get_kd(kd0, ΔHkd, temperature)
#         return kd0 * exp(-ΔHkd / (PolymerMembranes.R_J_MOL_K * temperature))
#     end
#     function get_kd(M, Nkd, ΔHkd, temperature, critical_temperature)
#         kd0 = get_kd0(M, Nkd, critical_temperature)
#         kd = get_kd(kd0, ΔHkd, temperature)
#         return kd
#     end
#     function get_b(M, Nb, ΔHb, temperature, critical_temperature)
#         b0 = get_b0(M, Nb, critical_temperature)
#         b = get_b(b0, ΔHb, temperature)
#         return b
#     end

#     function recover_ch(isotherm, b, kd, use_fugacity)
#         function sub_target_for_ch(singleton_ch_vec)
#             return rss(DualModeModel(singleton_ch_vec[1], b, kd), isotherm; use_fugacity=use_fugacity)
#         end
#         return Optim.minimizer(Optim.optimize(sub_target_for_ch, [5.], LBFGS()))[1]
#     end

# end

# function rss(constraints::VantHoffDualModeModelConstraints, isotherms::AbstractVector{<:IsothermData}, use_fugacity)
#     temps = [temperature(isotherm) for isotherm in isotherms]
#     bs = [VHDMHelper.get_b(constraints.b0, constraints.ΔHb, temp) for temp in temps] 
#     kds = [VHDMHelper.get_kd(constraints.kd0, constraints.ΔHkd, temp) for temp in temps]
#     chs = MVector{length(isotherms), Float64}(fill(0, length(isotherms)))

#     for idx in eachindex(bs, kds, isotherms)
#         b = bs[idx]
#         kd = kds[idx]
#         function sub_target_for_ch(singleton_ch_vec)
#             return rss(DualModeModel(singleton_ch_vec[1], b, kd), isotherms[idx]; use_fugacity=use_fugacity)
#         end
#         chs[idx] = Optim.minimizer(Optim.optimize(sub_target_for_ch, [5.], LBFGS()))[1]
#     end
    
#     dm_models = [DualModeModel(chs[idx], bs[idx], kds[idx]) for idx in eachindex(chs, bs, kds)]
#     dual_model_error = sum([rss(dm_models[idx], isotherms[idx]; use_fugacity=use_fugacity) for idx in eachindex(dm_models, isotherms)])
#     return dual_model_error
# end

# function rss(constraints::MultiGasVantHoffDualModeModelConstraints, isotherm_vectors::AbstractVector{<:AbstractVector{<:IsothermData}}, critical_temperatures::AbstractVector{<:Number}, use_fugacity)
    
#     total_rss = 0
#     for idx in eachindex(isotherm_vectors, critical_temperatures, constraints.ΔHbs, constraints.ΔHkds)
#         b0 = VHDMHelper.get_b0(constraints.M, constraints.Nb, critical_temperatures[idx]) 
#         kd0 = VHDMHelper.get_kd0(constraints.M, constraints.Nkd, critical_temperatures[idx])
#         individual_gas_constraints = VantHoffDualModeModelConstraints(constraints.ΔHbs[idx], constraints.ΔHkds[idx], b0, kd0)
#         total_rss += rss(individual_gas_constraints, isotherm_vectors[idx], use_fugacity)
#     end
#     return total_rss
# end


# function fit_multiple_gasses_and_temperatures_to_dualmode(critical_temps, corresponding_isotherms::AbstractVector{<:AbstractVector{<:IsothermData}}; use_fugacity=false) 
#     num_gasses = length(corresponding_isotherms)
#     if num_gasses != length(critical_temps)
#         throw(AssertionError("Number of critical temperatures did not match the top level of the corresponding isotherms."))
#     end
#     # isotherms_for_each_gas = [length(corresponding_isotherms[idx]) for idx in eachindex(corresponding_isotherms)]
#     # ch_vector_start_positions = append!([1], [sum(isotherms_for_each_gas[1:idx]) + 1 for idx in eachindex(isotherms_for_each_gas)[1:end-1]])
#     # ch_vector_end_positions = [sum(isotherms_for_each_gas[1:idx]) for idx in eachindex(isotherms_for_each_gas)]
    
#     function target(param_vec)
#         M = param_vec[1]
#         Nb = param_vec[2]
#         Nkd = param_vec[3]
#         position = 4
#         hb_values = param_vec[position:position + num_gasses - 1]
#         position += num_gasses
#         hkd_values = param_vec[position:position + num_gasses - 1]
#         # position += num_gasses
#         # ch_values = param_vec[position:end]
#         constraints = MultiGasVantHoffDualModeModelConstraints(M, Nb, Nkd, hb_values, hkd_values)

#         # ch_values_vectors = [ch_values[ch_vector_start_positions[idx]:ch_vector_end_positions[idx]] for idx in eachindex(ch_vector_start_positions, ch_vector_end_positions)]
#         error = rss(constraints, corresponding_isotherms, critical_temps, use_fugacity)
#         if isnan(error)
#             @show param_vec
#         end
#         return error
#     end

#     # initial_ch_value_vectors = []
#     initial_hb_values = []
#     initial_hkd_values = []
    
#     # fit the initial models and use them as the initial parameters for the general optimization
#     for isotherm_vector in corresponding_isotherms
#         standalone_models_for_gas = [fit_model(DualMode(), isotherm; use_fugacity=use_fugacity) for isotherm in isotherm_vector]
        
#         # chs = [initial_model.ch for initial_model in standalone_models_for_gas] 
#         # append!(initial_ch_value_vectors, chs)

#         bs = [initial_model.b for initial_model in standalone_models_for_gas] 
#         kds = [initial_model.kd for initial_model in standalone_models_for_gas]

#         temps = [temperature(isotherm) for isotherm in isotherm_vector]
#         inverse_temps = 1 ./ temps
#         log_bs = log.(bs)
#         log_kds = log.(kds)
#         negative_ΔHb_per_R, log_b0 = fit_linear_data(inverse_temps, log_bs)
#         negative_ΔHkd_per_R, log_kd0 = fit_linear_data(inverse_temps, log_kds)
#         ΔHb = -negative_ΔHb_per_R * R_J_MOL_K
#         ΔHkd = -negative_ΔHkd_per_R * R_J_MOL_K
#         append!(initial_hb_values, ΔHb)
#         append!(initial_hkd_values, ΔHkd)
#     end

#     initial_hb_values = min.(-1000, strip_measurement_to_value(initial_hb_values))
#     initial_hkd_values = min.(-1000, strip_measurement_to_value(initial_hkd_values))
#     # initial_ch_value_vectors = max.(1e-5, strip_measurement_to_value(initial_ch_value_vectors))
#     # @show initial_hb_values 
#     # @show initial_hkd_values 
#     # @show initial_ch_value_vectors 
#     # follows M, Nb, Nkd, init_hbs, init_hkds, 
#     initial_param_vec = strip_measurement_to_value([0.5, 1, 1, initial_hb_values..., initial_hkd_values...])
#     upper_bound = strip_measurement_to_value([1., 50., 50., fill(0, size(initial_hb_values))..., fill(0, size(initial_hkd_values))...])
#     lower_bound = strip_measurement_to_value([0., -50., -50., fill(-100000, size(initial_hb_values))..., fill(-100000, size(initial_hkd_values))...])
#     # res = bboptimize(target; SearchRange = collect(zip(lower_bound, upper_bound)))
#     # res = HiddenVHDMOptimizer.optimize_function(target, lower_bound, upper_bound)
#     # println(res)
#     # param_vec = best_fitness(res)

#     # res = Optim.optimize(target, lower_bound, upper_bound, initial_param_vec, Fminbox(BFGS()))
#     # res = Optim.optimize(target, initial_param_vec, ParticleSwarm(;lower = lower_bound, upper = upper_bound, n_particles = 5))
    
#     # res = Optim.optimize(target, initial_param_vec, NewtonTrustRegion())
#     println("Beginning SAMIN routine with initial: " * string(initial_param_vec))
#     res = Optim.optimize(target, lower_bound, upper_bound, initial_param_vec, SAMIN(; rt=0.6), 
#         Optim.Options(iterations = 10^7))
#     # println(res)
#     param_vec = Optim.minimizer(res)

#     # @show param_vec
# end


