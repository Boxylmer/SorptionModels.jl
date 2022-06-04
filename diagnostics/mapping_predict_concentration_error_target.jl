using Revise
using SorptionModels
using MembraneEOS
using MembraneBase
using Plots


char_co2 = [630, 300, 1.515, 44]
char_ch4 = [250, 215, 0.500, 16.04]
char_n2 = [160, 145, 0.943, 28.01]

tpbo_25_density = 1.393  # g/cm3

co2_bulk_phase = SL(char_co2...)
char_tpbo_valerio = [474, 900, 1.6624, 100000]
co2_tpbo25_phase_valerio = SL([474, 630], [900, 300], [1.6624, 1.515], [100000, 44], [0 0; 0 0])
co2_tpbo25_nelf_valerio = NELFModel(co2_bulk_phase, co2_tpbo25_phase_valerio, tpbo_25_density)
target = SorptionModels._make_nelf_model_mass_fraction_target(co2_tpbo25_nelf_valerio, 308.15, 1, [1])



# high res, one image
pstars_highres = 400:10:1300
tstars_highres = 400:10:1300
# rhostars_highres = [char_tpbo25[3]]
rhostars_highres = [1.9]
function make_char_param_error_map(target_func, added_text = "")
    needed_iters_highres = length(pstars_highres) * length(tstars_highres) * length(rhostars_highres)
    errs_2D = Array{Float64, 2}(undef, length(pstars_highres), length(tstars_highres))
    for p_i in eachindex(pstars_highres)
        pstar = pstars_highres[p_i]
        work_done = "Single image " * string(round(length(tstars_highres) * length(rhostars_highres) * p_i / needed_iters_highres * 100))*"% complete."
        println(work_done)
        Threads.@threads for t_i in eachindex(tstars_highres)
            tstar = tstars_highres[t_i]
            rhostar = rhostars_highres[1]
            errs_2D[p_i, t_i] = target_func([pstar, tstar, rhostar, mw])
        end
    end
    # err_2d_heatmap = heatmap(tstars_highres, pstars_highres, errs_2D[:, :, 1], title="rho="*string(rhostars_highres[1])*"g/cm3", xlabel = "T* (K)", ylabel = "P* (MPa)")
    err_2d_heatmap = contourf(
        tstars_highres, pstars_highres, errs_2D[:, :, 1], 
        title="rho="*string(rhostars_highres[1])*"g/cm3" * " " * added_text, 
        xlabel = "T* (K)", ylabel = "P* (MPa)")
    return err_2d_heatmap
end
err_2d_heatmap = make_char_param_error_map(error_target, "original")
err_2d_heatmap_2 = make_char_param_error_map(error_target_2, "new")
combined_plot = plot(err_2d_heatmap, err_2d_heatmap_2, layout = (1, 2))
savefig(combined_plot, joinpath(@__DIR__, "2D error heatmap (high res).png")) 


# # low res, animation
# pstars_lowres = 500:20:1300
# tstars_lowres = 500:20:1300 
# rhostars_lowres = 1.4:0.07:2.4

# function make_char_param_error_map_3d(target_func)
#     needed_iters_lowres = length(pstars_lowres) * length(tstars_lowres) * length(rhostars_lowres)
#     errs_3D = Array{Float64, 3}(undef, length(pstars_lowres), length(tstars_lowres), length(rhostars_lowres))
#     for p_i in eachindex(pstars_lowres)
#         pstar = pstars_lowres[p_i]
#         work_done = "Animation " * string(round(length(tstars_lowres) * length(rhostars_lowres) * p_i / needed_iters_lowres * 100))*"% complete."
#         println(work_done)
#         Threads.@threads for t_i in eachindex(tstars_lowres)
#             tstar = tstars_lowres[t_i]
#             for r_i in eachindex(rhostars_lowres)
#                 rhostar = rhostars_lowres[r_i]
#                 errs_3D[p_i, t_i, r_i] = target_func([pstar, tstar, rhostar, mw])
#             end
#         end
#     end
#     anim = @animate for i ∈ 1:size(errs_3D)[3]
#         contourf(tstars_lowres, pstars_lowres, errs_3D[:, :, i], title="rho="*string(rhostars_lowres[i])*"g/cm3",clim=(minimum(errs_3D),maximum(errs_3D)), xlabel = "T* (K)", ylabel = "P* (MPa)")
#     end
#     return anim
# end
# gif(make_char_param_error_map_3d(error_target_2), joinpath(@__DIR__, "3D animation of SL char vals.gif"), fps = 5)