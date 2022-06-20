using Revise
using SorptionModels
using MembraneBase
using Plots


partial_pres = [0.00011198, 0.000247609, 0.000521334, 0.000821977, 0.00114812]
acts = [0.035361112, 0.078190268, 0.164627403, 0.25956501, 0.362554938]
concs = [0.023724367, 0.040262462, 0.060287035, 0.072594922, 0.079677023]
iso_low_conc = IsothermData(partial_pressures_mpa=partial_pres, concentrations_cc=concs, activities=acts)
model_low_conc = fit_model(GAB(), iso_low_conc)
conc_pred = predict_concentration(model_low_conc, partial_pres)

error_target = SorptionModels._make_gab_target(acts, concs)

# high res, one image
cp_highres = 200:30:1200
k_highres = 200:30:1200
a_highres = [4233.4]
# rhostars_highres = [1.6624]
function make_error_map(target_func, added_text = "")
    needed_iters_highres = length(cp_highres) * length(k_highres) * length(a_highres)
    errs_2D = Array{Float64, 2}(undef, length(cp_highres), length(k_highres))
    for cp_i in eachindex(cp_highres)
        cp_val = cp_highres[cp_i]
        work_done = "Single image " * string(round(length(k_highres) * length(a_highres) * cp_i / needed_iters_highres * 100))*"% complete."
        println(work_done)
        # Threads.@threads 
        Threads.@threads for k_i in eachindex(k_highres)
            k_val = k_highres[k_i]
            a_val = a_highres[1]
            errs_2D[cp_i, k_i] = target_func([cp_val, k_val, a_val])
        end
    end
    # err_2d_heatmap = heatmap(tstars_highres, pstars_highres, errs_2D[:, :, 1], title="rho="*string(rhostars_highres[1])*"g/cm3", xlabel = "T* (K)", ylabel = "P* (MPa)")
    err_2d_heatmap = contourf(
        k_highres, cp_highres, errs_2D[:, :, 1], 
        title="rho="*string(round(a_highres[1]; digits=2))*"g/cm3" * " " * added_text, 
        xlabel = "T* (K)", ylabel = "P* (MPa)", size=(600,600), clims=(minimum(errs_2D), sum(errs_2D)/length(errs_2D)))
    return err_2d_heatmap
end
err_2d_heatmap = make_error_map(error_target, "original", )
# err_2d_heatmap_2 = make_char_param_error_map(error_target_2, "sinf + isotherm")
# err_2d_heatmap_3 = make_char_param_error_map(error_target_3, "sarti 2001")
# combined_plot = plot(err_2d_heatmap, err_2d_heatmap_3, layout = 2)
# savefig(combined_plot, joinpath(@__DIR__, "2D error heatmap (high res).png"))


# low res, animation
cp_lowres = 0:20:1400
k_lowres = 0:20:1400
a_lowres = 0:20:100

function make_char_param_error_map_3d(target_func)
    needed_iters_lowres = length(cp_lowres) * length(k_lowres) * length(a_lowres)
    errs_3D = Array{Float64, 3}(undef, length(cp_lowres), length(k_lowres), length(a_lowres))
    for cp_i in eachindex(cp_lowres)
        cp_val = cp_lowres[cp_i]
        work_done = "Animation " * string(round(length(k_lowres) * length(a_lowres) * cp_i / needed_iters_lowres * 100))*"% complete."
        println(work_done)
        Threads.@threads for k_i in eachindex(k_lowres)
            k_val = k_lowres[k_i]
            for a_i in eachindex(a_lowres)
                a_val = a_lowres[a_i]
                errs_3D[cp_i, k_i, a_i] = target_func([cp_val, k_val, a_val, mw])
            end
        end
    end

    minlocation = argmin(errs_3D)
    @show global_min = [k_lowres[minlocation[2]]], [cp_lowres[minlocation[1]]]
    frameminlocations = [argmin(errs_3D[:, :, i]) for i in eachindex(a_lowres)]
    mintstarvals = [k_lowres[frameminlocations[i][2]] for i in eachindex(a_lowres)]
    minpstarvals = [cp_lowres[frameminlocations[i][1]] for i in eachindex(a_lowres)]

    anim = @animate for i ∈ 1:size(errs_3D)[3]
        myplot = contourf(k_lowres, cp_lowres, errs_3D[:, :, i], title="rho="*string(round(a_lowres[i]; digits=3))*"g/cm3",clim=(minimum(errs_3D),5), xlabel = "T* (K)", ylabel = "P* (MPa)")
        annotate!(myplot, (0.8, 0.9), "min =" * string(round(minimum(errs_3D[:, :, i]); digits=3)))
        scatter!(myplot, [mintstarvals], [minpstarvals], legend=false)
        scatter!(myplot, global_min..., legend=false, markersize=6)
    end
    return anim
end
# gif(make_char_param_error_map_3d(error_target), joinpath(@__DIR__, "3D animation of SL char vals.gif"), fps = 2)
# gif(make_char_param_error_map_3d(error_target_3), joinpath(@__DIR__, "3D animation of SL char vals, sarti plot.gif"), fps = 5)