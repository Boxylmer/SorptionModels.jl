using Optim
using DifferentialEquations
using Plots
using LinearAlgebra
using MembraneBase


# Define the derivative function
struik_dfdt(f, fe, τ∞, γ) = -(f - fe) / (τ∞ * exp(-γ * (f - fe)))

function struik_dfdt!(df, f, p, _)
    fe, τ∞, γ = p
    df[1] = struik_dfdt(f[1], fe, τ∞, γ)
end

function solve_struik(fe, τ∞, γ, f0, times)
    p = (fe, τ∞, γ)
    prob = ODEProblem(struik_dfdt!, [f0], (0, times[end]), p)
    sol = solve(prob, saveat=times)
    return [v[1] for v in sol.u]
end

# test it
solve_struik(0.1, 0.1, 0.1, 1.0, collect(range(0, 10, 100)))

function struik_sse(times, fs, params)
    f0 = fs[1]
    pred = solve_struik(params..., f0, times)
    
    if length(pred) != length(fs)
        return NaN
    else
        return sum((pred .- fs) .^ 2)
    end

end

# function struik_sse_adjustable_f0(times, fs, params)
#     pred = solve_struik(params..., times)
#     if length(pred) != length(fs)
#         return NaN
#     else
#         return sum((pred .- fs) .^ 2)
#     end
# end

times = collect(range(0, 10, 10))
fs = [0.4, 0.3, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
struik_sse(times, fs, [0.2, 150, 130])
# struik_sse_adjustable_f0(times, fs, [0.2, 150, 130, 0.21])


function solve_struik_params(times, fs)
    function obj(params)
        lnfe, lnτ∞, lnγ = (params[1], params[2], params[3])
        err = log(
            struik_sse(
                times, 
                fs, 
                exp.((lnfe, lnτ∞, lnγ))
            )
        )
        return err
    end
    start = [log(fs[end]), log(100), log(10)]
    res = Optim.optimize(obj, start, NelderMead()).minimizer
    lnfe, lnτ∞, lnγ = res
    params = exp.((lnfe, lnτ∞, lnγ))
    return params
end


# function solve_struik_params_adjusable_f0(times, fs)
#     function obj(params)
#         fe, lnτ∞, γ, f0 = (params[1], params[2], params[3], params[4])
#         err = log(struik_sse_adjustable_f0(times, fs, (fe , exp(lnτ∞), γ, f0)))
#         return err
#     end
#     start = [fs[end], log(100), log(100), fs[1]]
#     res = Optim.optimize(obj, start, NelderMead()).minimizer
#     fe, lnτ∞, γ, f0 = res
#     params = (fe, exp(lnτ∞), γ, f0)
#     return params
# end


# real data times
f_from_p(p, A = 57663.53, B = 0.969509) = -B / log(p/A)
# NEAT_PTMSP_TIMES_H = [0.677, 0.828, 0.914, 0.969, 1.685, 1.908, 4.831, 5.000, 5.665, 5.925, 6.692361111, 6.88125, 7.752777778, 8.877083333, 13.95, 14.88263889, 15.92916667, 19.68194444, 20.752, 21.659] .* 24
NEAT_PTMSP_TIMES_H = [0.000, 3.633, 5.683, 7.000, 24.183, 28.750, 99.700, 103.750, 119.717, 125.950, 144.367, 148.900, 169.817, 196.800, 318.550, 340.933, 366.050, 456.117, 481.800, 503.567]
NEAT_PTMSP_PERMS_BARRER = [822.153, 810.529, 799.437, 798.998, 753.495, 745.478, 658.444, 654.332, 643.441, 638.246, 623.815, 622.938, 609.143, 598.718, 544.091, 535.868, 524.821, 500.271, 489.880, 487.181]
NEAT_PTMSP_FS = f_from_p.(NEAT_PTMSP_PERMS_BARRER)

NEAT_PTMSP_PARAMS = solve_struik_params(NEAT_PTMSP_TIMES_H, NEAT_PTMSP_FS)
NEAT_PTMSP_PRED = solve_struik(NEAT_PTMSP_PARAMS..., NEAT_PTMSP_FS[1], NEAT_PTMSP_TIMES_H)
# NEAT_PTMSP_PARAMS_WITH_F0 = solve_struik_params_adjusable_f0(NEAT_PTMSP_TIMES_H, NEAT_PTMSP_FS)
# NEAT_PTMSP_PRED_WITH_F0 = solve_struik(NEAT_PTMSP_PARAMS_WITH_F0..., NEAT_PTMSP_TIMES_H)


# see the error map of struik_res with gamma and tau inf
function generate_τγ_error_map(fe=0.19, lnτ∞_range = 5:10:100, γ_range = 1:0.05:2)
    
    error_matrix = zeros(length(γ_range), length(lnτ∞_range))
    
    for i in eachindex(γ_range)
        for j in eachindex(lnτ∞_range)
            error_matrix[i, j] = log(struik_sse(
                NEAT_PTMSP_TIMES_H, 
                NEAT_PTMSP_FS, 
                (fe, exp(lnτ∞_range[j]), γ_range[i])))
        end
    end
    
    heatmap(γ_range, lnτ∞_range, error_matrix',
            xlabel="γ",
            ylabel="lnτ∞",
            title="Error Heatmap fe=" * string(round(fe, digits=2)),
            color=:viridis)  # Use a color scheme that suits your data; :viridis is a good default
end  

function generate_feγ_error_map(fe_range=0.1:0.001:0.2, τ∞=5, γ_range=1:0.05:2)
    error_matrix = zeros(length(γ_range), length(fe_range))
    
    for i in eachindex(γ_range)
        for j in eachindex(fe_range)
            error_matrix[i, j] = log(
                struik_sse(
                    NEAT_PTMSP_TIMES_H, 
                    NEAT_PTMSP_FS, 
                    (fe_range[j], τ∞, γ_range[i])
                )
            )
        end
    end
    
    heatmap(γ_range, fe_range, error_matrix',
            xlabel="γ",
            ylabel="fe",
            title="Error Heatmap τ∞=" * string(round(τ∞, digits=2)),
            color=:viridis)  # Use a color scheme that suits your data; :viridis is a good default
end  

function generate_feτ_error_map(fe_range=0.1:0.001:0.2, lnτ∞_range=5:10:100, γ=10)
    error_matrix = zeros(length(fe_range), length(lnτ∞_range))
    @show length(lnτ∞_range), length(fe_range)
    for i in eachindex(fe_range)
        for j in eachindex(lnτ∞_range)
            error_matrix[i, j] = log(
                struik_sse(
                    NEAT_PTMSP_TIMES_H, 
                    NEAT_PTMSP_FS, 
                    (fe_range[i], exp(lnτ∞_range[j]), γ)
                )
            )
        end
    end
    
    heatmap(fe_range, lnτ∞_range, error_matrix',
            xlabel="fe",
            ylabel="lnτ∞",
            title="Error Heatmap γ=" * string(round(γ, digits=2)),
            color=:viridis)  # Use a color scheme that suits your data; :viridis is a good default
end  

function generate_feτγ_error_map(fe_range, lnτ∞_range, γ_range, given_times, given_fs)
    
    error_matrix = zeros(length(γ_range), length(lnτ∞_range), length(fe_range))
    
    for i in eachindex(γ_range)
        for j in eachindex(lnτ∞_range)
            for k in eachindex(fe_range)

                error_matrix[i, j, k] = log(struik_sse(
                    given_times, 
                    given_fs, 
                    (fe, exp(lnτ∞_range[j]), γ_range[i])))
            end
        end
    end
    
    heatmap(γ_range, lnτ∞_range, error_matrix',
            xlabel="γ",
            ylabel="lnτ∞",
            title="Error Heatmap fe=" * string(round(fe, digits=2)),
            color=:viridis)  # Use a color scheme that suits your data; :viridis is a good default
end  

τγ_plot = generate_τγ_error_map(NEAT_PTMSP_PARAMS[1], 8:0.1:32, 10:1:200)  
feγ_plot = generate_feγ_error_map(0.01:0.001:0.20, NEAT_PTMSP_PARAMS[2], 10:1:200)
feτ_plot = generate_feτ_error_map(0.01:0.001:0.20, 8:0.1:32, NEAT_PTMSP_PARAMS[3])
cumulative_plot = plot(τγ_plot, feγ_plot, feτ_plot, layout=[2, 2], legend=false)


function generate_fτγ_error_animation(fe_range=0.09:0.01:0.22, lnτ∞_range = 5:10:100, γ_range = 5:10:100; name = "fτγ_error_animation.gif")
    anim = @animate for (i, fe) in enumerate(fe_range)
        print("\r" * string(round(i/length(fe_range), digits=3)))
        generate_τγ_error_map(fe, lnτ∞_range, γ_range)
    end
    return gif(anim, name, fps = 8)
end 
generate_fτγ_error_animation(0.05:0.01:0.21, 8:0.2:32, 0:1:200) 

plt = scatter(NEAT_PTMSP_TIMES_H, NEAT_PTMSP_FS)
prediction_times = collect(range(0, NEAT_PTMSP_TIMES_H[end], 250))
plt = plot!(prediction_times, solve_struik(NEAT_PTMSP_PARAMS..., NEAT_PTMSP_FS[1], prediction_times))
# plt = plot!(prediction_times, solve_struik(NEAT_PTMSP_PARAMS_WITH_F0..., prediction_times))

NEAT_PTMSP_PARAMS_Σ = MembraneBase.rss_covariance_matrix(
    params -> struik_sse(NEAT_PTMSP_TIMES_H, NEAT_PTMSP_FS, params), 
    collect(NEAT_PTMSP_PARAMS), 
    length(NEAT_PTMSP_FS)
) 
NEAT_PTMSP_PARAMS_σ = sqrt.(abs.(diag(NEAT_PTMSP_PARAMS_Σ)))
generate_τγ_error_map(NEAT_PTMSP_PARAMS[1], 8:0.1:32, 10:1:200)  


