using Revise
using SorptionModels
using MembraneEOS
using MembraneBase
using Plots

# test the polymer fitter with TPBO-0.25
tpbo_ch4_5c = IsothermData(; 
    partial_pressures_mpa = [0, 0.051022404, 0.097858902, 0.172285293, 0.272711361, 0.371012386, 0.475068139], 
    concentrations_cc = [0, 7.37368523, 12.0433614, 17.76552202, 23.28373709, 27.50367509, 31.07457011],
    temperature_k = 278.15, rho_pol_g_cm3 = 1.3937
    )
tpbo_ch4_35c = IsothermData(; 
    partial_pressures_mpa = [0, 0.05676287, 0.103720596, 0.177868877, 0.268361442, 0.371119351, 0.478013248], 
    concentrations_cc = [0, 3.741224553, 6.311976164, 9.748565324, 13.23714075, 16.47955269, 19.49981169],
    temperature_k = 308.15, rho_pol_g_cm3 = 1.3937
    )
tpbo_co2_20c = IsothermData(; 
    partial_pressures_mpa = [0, 0.018051448, 0.043568055, 0.07237145, 0.134854673, 0.239969739, 0.386544681], 
    concentrations_cc = [0, 11.64079279, 22.18344653, 30.5524273, 43.09494749, 56.42684262, 67.28267947],
    temperature_k = 293.15, rho_pol_g_cm3 = 1.3937
    )
tpbo_co2_50c = IsothermData(; 
    partial_pressures_mpa = [0, 0.023513454, 0.050773712, 0.080001807, 0.144376557, 0.249710838, 0.396483131], 
    concentrations_cc = [0, 5.93630284, 11.36554572, 15.98552528, 23.62447856, 33.06426088, 42.47173453],
    temperature_k = 323.15, rho_pol_g_cm3 = 1.3937
    )
tpbo_n2_5c = IsothermData(; 
    partial_pressures_mpa = [0, 0.068906986, 0.181377336, 0.424951374, 0.731306858, 1.064696014, 1.413103086], 
    concentrations_cc = [0, 2.252715738, 5.601581157, 11.31054253, 16.87930294, 21.39238669, 25.17075548],
    temperature_k = 278.15, rho_pol_g_cm3 = 1.3937
    )
tpbo_n2_50c = IsothermData(; 
    partial_pressures_mpa = [0, 0.269930265, 0.705742173, 1.060688385, 1.42192415, 1.813024602, 2.228349107], 
    concentrations_cc = [0, 2.435212223, 5.677614879, 8.139676474, 10.6450967, 12.90356804, 14.82380991],
    temperature_k = 323.15, rho_pol_g_cm3 = 1.3937
    )
isotherms = [tpbo_ch4_5c, tpbo_ch4_35c, tpbo_co2_20c, tpbo_co2_50c, tpbo_n2_5c, tpbo_n2_50c]

char_co2 = [630, 300, 1.515, 44]
char_ch4 = [250, 215, 0.500, 16.04]
char_n2 = [160, 145, 0.943, 28.01]

tpbo_25_density = 1.393  # g/cm3

tpbo_co2_50c_exp = concentration(tpbo_co2_50c; component=1)
tpbo_co2_20c_exp = concentration(tpbo_co2_20c; component=1)
co2_bulk_phase = SL(char_co2...)
bulk_phase_char_params = [char_ch4, char_ch4, char_co2, char_co2, char_n2, char_n2]

char_tpbo_valerio = [474, 900, 1.6624, 100000]
co2_tpbo25_phase_valerio = SL([474, 630], [900, 300], [1.6624, 1.515], [100000, 44], [0 0; 0 0])
co2_tpbo25_nelf_valerio = NELFModel(co2_bulk_phase, co2_tpbo25_phase_valerio, tpbo_25_density)
tpbo_co2_50c_valerio = [predict_concentration(co2_tpbo25_nelf_valerio, 323.15, p)[1] for p in partial_pressures(tpbo_co2_50c; component=1)]
tpbo_co2_20c_valerio = [predict_concentration(co2_tpbo25_nelf_valerio, 293.15, p)[1] for p in partial_pressures(tpbo_co2_20c; component=1)]

#fit char params
char_tpbo25 = fit_model(NELF(), isotherms, bulk_phase_char_params)
@show char_tpbo25
co2_tpbo25_phase_fitted = SL([char_tpbo25[1], 630], [char_tpbo25[2], 300], [char_tpbo25[3], 1.515], [char_tpbo25[4], 44], [0 0.0; 0.0 0])
co2_tpbo25_nelf_fitted = NELFModel(co2_bulk_phase, co2_tpbo25_phase_fitted, tpbo_25_density)
tpbo_co2_50c_fitted = [predict_concentration(co2_tpbo25_nelf_fitted, 323.15, p)[1] for p in partial_pressures(tpbo_co2_50c; component=1)]
tpbo_co2_20c_fitted = [predict_concentration(co2_tpbo25_nelf_fitted, 293.15, p)[1] for p in partial_pressures(tpbo_co2_20c; component=1)]

tpbo_co2_50c_plot = plot(partial_pressures(tpbo_co2_50c; component=1), tpbo_co2_50c_valerio, label="valerio 50C", legend=:topleft)
plot!(tpbo_co2_50c_plot, partial_pressures(tpbo_co2_50c; component=1), tpbo_co2_50c_fitted, label="fitted 50C")
plot!(tpbo_co2_50c_plot, partial_pressures(tpbo_co2_50c; component=1), tpbo_co2_50c_exp, label="exp 50C")
tpbo_co2_20c_plot = plot(partial_pressures(tpbo_co2_20c; component=1), tpbo_co2_20c_valerio, label="valerio 20C", legend=:topleft)
plot!(tpbo_co2_20c_plot, partial_pressures(tpbo_co2_20c; component=1), tpbo_co2_20c_fitted, label="fitted 20C")
plot!(tpbo_co2_20c_plot, partial_pressures(tpbo_co2_20c; component=1), tpbo_co2_20c_exp, label="exp 20C")

tpbo_co2_plot = plot(tpbo_co2_50c_plot, tpbo_co2_20c_plot, layout=(2, 1))

savefig(tpbo_co2_plot, joinpath(@__DIR__, "tpbo_co2_50c_plot_comparison.png"))



error_target = SorptionModels._make_nelf_model_parameter_target(isotherms, bulk_phase_char_params, 1e-5)
error_target_2 = SorptionModels._make_nelf_model_parameter_target_2(isotherms, bulk_phase_char_params, 1e-5)
# @show error_target(char_tpbo25)
# @show error_target(char_tpbo_valerio)
return 3
mw = 1e6

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


# low res, animation
pstars_lowres = 500:20:1300
tstars_lowres = 500:20:1300 
rhostars_lowres = 1.4:0.07:2.4

function make_char_param_error_map_3d(target_func)
    needed_iters_lowres = length(pstars_lowres) * length(tstars_lowres) * length(rhostars_lowres)
    errs_3D = Array{Float64, 3}(undef, length(pstars_lowres), length(tstars_lowres), length(rhostars_lowres))
    for p_i in eachindex(pstars_lowres)
        pstar = pstars_lowres[p_i]
        work_done = "Animation " * string(round(length(tstars_lowres) * length(rhostars_lowres) * p_i / needed_iters_lowres * 100))*"% complete."
        println(work_done)
        Threads.@threads for t_i in eachindex(tstars_lowres)
            tstar = tstars_lowres[t_i]
            for r_i in eachindex(rhostars_lowres)
                rhostar = rhostars_lowres[r_i]
                errs_3D[p_i, t_i, r_i] = target_func([pstar, tstar, rhostar, mw])
            end
        end
    end
    anim = @animate for i ∈ 1:size(errs_3D)[3]
        contourf(tstars_lowres, pstars_lowres, errs_3D[:, :, i], title="rho="*string(rhostars_lowres[i])*"g/cm3",clim=(minimum(errs_3D),maximum(errs_3D)), xlabel = "T* (K)", ylabel = "P* (MPa)")
    end
    return anim
end
gif(make_char_param_error_map_3d(error_target_2), joinpath(@__DIR__, "3D animation of SL char vals.gif"), fps = 5)