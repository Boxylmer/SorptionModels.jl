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
tpbo_ch4_20c = IsothermData(; 
    partial_pressures_mpa = [0, 0.054409635, 0.101117776, 0.17570818, 0.275518841, 0.373215869], 
    concentrations_cc = [0, 5.084334416, 8.57780684, 12.92845101, 17.55678063, 21.17207162],
    temperature_k = 293.15, rho_pol_g_cm3 = 1.3937
    )
tpbo_ch4_35c = IsothermData(; 
    partial_pressures_mpa = [0, 0.05676287, 0.103720596, 0.177868877, 0.268361442, 0.371119351, 0.478013248], 
    concentrations_cc = [0, 3.741224553, 6.311976164, 9.748565324, 13.23714075, 16.47955269, 19.49981169],
    temperature_k = 308.15, rho_pol_g_cm3 = 1.3937
    )
tpbo_co2_5c = IsothermData(; 
    partial_pressures_mpa = [0, 0.012474405, 0.043391439, 0.11294484, 0.260641173, 0.563818055, 1.023901167, 1.633212939, 2.105369806, 2.615898101], 
    concentrations_cc = [0, 15.19436518, 33.5283133, 54.43691918, 76.32730082, 99.6986546, 121.6846495, 142.8977897, 157.6456881, 174.186204],
    temperature_k = 278.15, rho_pol_g_cm3 = 1.3937
    )
tpbo_co2_20c = IsothermData(; 
    partial_pressures_mpa = [0, 0.018051448, 0.043568055, 0.07237145, 0.134854673, 0.239969739, 0.386544681], 
    concentrations_cc = [0, 11.64079279, 22.18344653, 30.5524273, 43.09494749, 56.42684262, 67.28267947],
    temperature_k = 293.15, rho_pol_g_cm3 = 1.3937
    )
tpbo_co2_35c = IsothermData(; 
    partial_pressures_mpa = [0, 0.031852099, 0.066294896, 0.104825903, 0.142384952, 0.18000105, 0.217614795, 0.29451087, 0.380689171], 
    concentrations_cc = [0, 11.7485186, 19.85668992, 26.48815715, 31.5333868, 36.0129803, 39.61060375, 45.73820821, 51.29191881],
    temperature_k = 308.15, rho_pol_g_cm3 = 1.3937
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
isotherms = [tpbo_ch4_5c, tpbo_ch4_20c, tpbo_ch4_35c, tpbo_co2_5c, tpbo_co2_20c, tpbo_co2_35c, tpbo_co2_50c, tpbo_n2_5c, tpbo_n2_50c]
# dualmode_models = [fit_model(DualMode(), isotherm) for isotherm in isotherms]


function plot_dualmode_sinf(isotherm::IsothermData, isotherm_dualmode_exp_comparison=plot())
    isotherm_dualmode = fit_model(DualMode(), isotherm)
    isotherm_pressures = partial_pressures(isotherm; component=1)
    scatter!(isotherm_dualmode_exp_comparison, isotherm_pressures, concentration(isotherm; component=1), )#label="exp")
    isotherm_max_p = maximum(isotherm_pressures)
    isotherm_pred_pressures = 0:0.005:isotherm_max_p
    isotherm_dualmode_pred = predict_concentration(isotherm_dualmode, isotherm_pred_pressures)
    plot!(isotherm_dualmode_exp_comparison, isotherm_pred_pressures, isotherm_dualmode_pred, )#label ="dualmode")

    isotherm_sinf_empirical_pressure = 1e-5
    isotherm_sinf = predict_concentration(isotherm_dualmode, isotherm_sinf_empirical_pressure) / isotherm_sinf_empirical_pressure
    isotherm_sinf = SorptionModels.infinite_dilution_solubility(isotherm_dualmode)

    isotherm_sinf_pressures = [0, isotherm_max_p / 2]
    isotherm_sinf_concs = (isotherm_sinf .* isotherm_sinf_pressures)
    plot!(isotherm_dualmode_exp_comparison, isotherm_sinf_pressures, isotherm_sinf_concs, ) # label="infinite dilution slope")
end
plot_dualmode_sinf!(myplot, isotherm) = plot_dualmode_sinf(isotherm, myplot)
function plot_dualmode_sinf(isotherms::AbstractVector)
    isotherm_dualmode_exp_comparison=plot()
    for isotherm in isotherms
        plot_dualmode_sinf!(isotherm_dualmode_exp_comparison, isotherm)
    end
    return isotherm_dualmode_exp_comparison
end
tpbo_co2_comparison_plot = plot_dualmode_sinf([ tpbo_co2_20c, tpbo_co2_35c, tpbo_co2_50c])
savefig(tpbo_co2_comparison_plot, joinpath(@__DIR__, "tpbo_co2_dualmode_sinf_verification.png"))

char_co2 = [630, 300, 1.515, 44]
char_ch4 = [250, 215, 0.500, 16.04]
char_n2 = [160, 145, 0.943, 28.01]

tpbo_25_density = 1.393  # g/cm3

tpbo_co2_50c_exp = concentration(tpbo_co2_50c; component=1)
tpbo_co2_35c_exp = concentration(tpbo_co2_35c; component=1)
tpbo_co2_20c_exp = concentration(tpbo_co2_20c; component=1)
tpbo_co2_5c_exp = concentration(tpbo_co2_5c; component=1)

co2_bulk_phase = SL(char_co2...)
bulk_phase_char_params = [char_ch4, char_ch4, char_ch4, char_co2, char_co2, char_co2, char_co2, char_n2, char_n2]

char_tpbo_valerio = [474, 900, 1.6624, 100000]
co2_tpbo25_phase_valerio = SL([474, 630], [900, 300], [1.6624, 1.515], [100000, 44], [0 0; 0 0])
co2_tpbo25_nelf_valerio = NELFModel(co2_bulk_phase, co2_tpbo25_phase_valerio, tpbo_25_density)
tpbo_co2_50c_valerio = [predict_concentration(co2_tpbo25_nelf_valerio, 323.15, p)[1] for p in partial_pressures(tpbo_co2_50c; component=1)]
tpbo_co2_35c_valerio = [predict_concentration(co2_tpbo25_nelf_valerio, 308.15, p)[1] for p in partial_pressures(tpbo_co2_35c; component=1)]
tpbo_co2_20c_valerio = [predict_concentration(co2_tpbo25_nelf_valerio, 293.15, p)[1] for p in partial_pressures(tpbo_co2_20c; component=1)]
tpbo_co2_5c_valerio = [predict_concentration(co2_tpbo25_nelf_valerio, 278.15, p)[1] for p in partial_pressures(tpbo_co2_5c; component=1)]

#fit char params
char_tpbo25 = fit_model(NELF(), isotherms, bulk_phase_char_params)
@show char_tpbo25
kij_test =  0.00
ksw_test = [0.00]
co2_tpbo25_phase_fitted = SL([char_tpbo25[1], 630], [char_tpbo25[2], 300], [char_tpbo25[3], 1.515], [char_tpbo25[4], 44], [0 kij_test; kij_test 0])
co2_tpbo25_nelf_fitted = NELFModel(co2_bulk_phase, co2_tpbo25_phase_fitted, tpbo_25_density)
tpbo_co2_50c_fitted = [predict_concentration(co2_tpbo25_nelf_fitted, 323.15, p; ksw=ksw_test)[1] for p in partial_pressures(tpbo_co2_50c; component=1)]
tpbo_co2_35c_fitted = [predict_concentration(co2_tpbo25_nelf_fitted, 308.15, p; ksw=ksw_test)[1] for p in partial_pressures(tpbo_co2_35c; component=1)]
tpbo_co2_20c_fitted = [predict_concentration(co2_tpbo25_nelf_fitted, 293.15, p; ksw=ksw_test)[1] for p in partial_pressures(tpbo_co2_20c; component=1)]
tpbo_co2_5c_fitted = [predict_concentration(co2_tpbo25_nelf_fitted, 278.15, p; ksw=ksw_test)[1] for p in partial_pressures(tpbo_co2_5c; component=1)]

tpbo_co2_50c_plot = plot(partial_pressures(tpbo_co2_50c; component=1), tpbo_co2_50c_valerio, label="valerio 50C", legend=:topleft)
plot!(tpbo_co2_50c_plot, partial_pressures(tpbo_co2_50c; component=1), tpbo_co2_50c_fitted, label="fitted 50C")
plot!(tpbo_co2_50c_plot, partial_pressures(tpbo_co2_50c; component=1), tpbo_co2_50c_exp, label="exp 50C")
tpbo_co2_35c_plot = plot(partial_pressures(tpbo_co2_35c; component=1), tpbo_co2_35c_valerio, label="valerio 35C", legend=:topleft)
plot!(tpbo_co2_35c_plot, partial_pressures(tpbo_co2_35c; component=1), tpbo_co2_35c_fitted, label="fitted 35C")
plot!(tpbo_co2_35c_plot, partial_pressures(tpbo_co2_35c; component=1), tpbo_co2_35c_exp, label="exp 35C")
tpbo_co2_20c_plot = plot(partial_pressures(tpbo_co2_20c; component=1), tpbo_co2_20c_valerio, label="valerio 20C", legend=:topleft)
plot!(tpbo_co2_20c_plot, partial_pressures(tpbo_co2_20c; component=1), tpbo_co2_20c_fitted, label="fitted 20C")
plot!(tpbo_co2_20c_plot, partial_pressures(tpbo_co2_20c; component=1), tpbo_co2_20c_exp, label="exp 20C")
tpbo_co2_5c_plot = plot(partial_pressures(tpbo_co2_5c; component=1), tpbo_co2_5c_valerio, label="valerio 5C", legend=:topleft)
plot!(tpbo_co2_5c_plot, partial_pressures(tpbo_co2_5c; component=1), tpbo_co2_5c_fitted, label="fitted 5C")
plot!(tpbo_co2_5c_plot, partial_pressures(tpbo_co2_5c; component=1), tpbo_co2_5c_exp, label="exp 5C")
tpbo_co2_plot = plot(tpbo_co2_50c_plot, tpbo_co2_35c_plot, tpbo_co2_20c_plot, tpbo_co2_5c_plot, layout=(2, 2))
# plot!(tpbo_co2_plot, xlims=(0, 0.1))
savefig(tpbo_co2_plot, joinpath(@__DIR__, "tpbo_co2_50c_plot_comparison.png"))

inf_dil_pres = 1e-10
error_target = SorptionModels._make_nelf_model_parameter_target(isotherms, bulk_phase_char_params, inf_dil_pres)
error_target_2 = SorptionModels._make_nelf_model_parameter_target_2(isotherms, bulk_phase_char_params, inf_dil_pres)
error_target_3 = SorptionModels._make_nelf_model_parameter_target_3(isotherms, bulk_phase_char_params, inf_dil_pres)
# @show error_target(char_tpbo25)
# @show error_target(char_tpbo_valerio)
mw = 1e6

# high res, one image
pstars_highres = 200:50:1200
tstars_highres = 200:50:1200
rhostars_highres = [char_tpbo25[3]]
# rhostars_highres = [1.6624]
function make_char_param_error_map(target_func, added_text = "")
    needed_iters_highres = length(pstars_highres) * length(tstars_highres) * length(rhostars_highres)
    errs_2D = Array{Float64, 2}(undef, length(pstars_highres), length(tstars_highres))
    for p_i in eachindex(pstars_highres)
        pstar = pstars_highres[p_i]
        work_done = "Single image " * string(round(length(tstars_highres) * length(rhostars_highres) * p_i / needed_iters_highres * 100))*"% complete."
        println(work_done)
        # Threads.@threads 
        Threads.@threads for t_i in eachindex(tstars_highres)
            tstar = tstars_highres[t_i]
            rhostar = rhostars_highres[1]
            errs_2D[p_i, t_i] = target_func([pstar, tstar, rhostar, mw])
        end
    end
    # err_2d_heatmap = heatmap(tstars_highres, pstars_highres, errs_2D[:, :, 1], title="rho="*string(rhostars_highres[1])*"g/cm3", xlabel = "T* (K)", ylabel = "P* (MPa)")
    err_2d_heatmap = contourf(
        tstars_highres, pstars_highres, errs_2D[:, :, 1], 
        title="rho="*string(round(rhostars_highres[1]; digits=2))*"g/cm3" * " " * added_text, 
        xlabel = "T* (K)", ylabel = "P* (MPa)", size=(600,600))
    return err_2d_heatmap
end
err_2d_heatmap = make_char_param_error_map(error_target, "sinf", )
err_2d_heatmap_2 = make_char_param_error_map(error_target_2, "sinf + isotherm")
err_2d_heatmap_3 = make_char_param_error_map(error_target_3, "sarti 2001")
combined_plot = plot(err_2d_heatmap, err_2d_heatmap_2, err_2d_heatmap_3, layout = 3,)
savefig(combined_plot, joinpath(@__DIR__, "2D error heatmap (high res).png"))


# low res, animation
pstars_lowres = 500:70:1200
tstars_lowres = 500:70:1200
rhostars_lowres = tpbo_25_density + 0.1:0.2:char_tpbo25[3] + 1

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
gif(make_char_param_error_map_3d(error_target), joinpath(@__DIR__, "3D animation of SL char vals.gif"), fps = 5)
# gif(make_char_param_error_map_3d(error_target_3), joinpath(@__DIR__, "3D animation of SL char vals, sarti plot.gif"), fps = 5)