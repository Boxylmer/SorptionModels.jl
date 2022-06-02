using SorptionModels
using MembraneEOS
using MembraneBase
using Plots
using Revise


# test the polymer fitter with TPBO-0.25
tpbo_ch4_5c = IsothermData(; 
    partial_pressures_mpa = [0.051022404, 0.097858902, 0.172285293, 0.272711361, 0.371012386, 0.475068139], 
    concentrations_cc = [7.37368523, 12.0433614, 17.76552202, 23.28373709, 27.50367509, 31.07457011],
    temperature_k = 278.15, rho_pol_g_cm3 = 1.3937
    )
tpbo_ch4_35c = IsothermData(; 
    partial_pressures_mpa = [0.05676287, 0.103720596, 0.177868877, 0.268361442, 0.371119351, 0.478013248], 
    concentrations_cc = [3.741224553, 6.311976164, 9.748565324, 13.23714075, 16.47955269, 19.49981169],
    temperature_k = 308.15, rho_pol_g_cm3 = 1.3937
    )
tpbo_co2_20c = IsothermData(; 
    partial_pressures_mpa = [0.018051448, 0.043568055, 0.07237145, 0.134854673, 0.239969739, 0.386544681], 
    concentrations_cc = [11.64079279, 22.18344653, 30.5524273, 43.09494749, 56.42684262, 67.28267947],
    temperature_k = 293.15, rho_pol_g_cm3 = 1.3937
    )
tpbo_co2_50c = IsothermData(; 
    partial_pressures_mpa = [0.023513454, 0.050773712, 0.080001807, 0.144376557, 0.249710838, 0.396483131], 
    concentrations_cc = [5.93630284, 11.36554572, 15.98552528, 23.62447856, 33.06426088, 42.47173453],
    temperature_k = 323.15, rho_pol_g_cm3 = 1.3937
    )
tpbo_n2_5c = IsothermData(; 
    partial_pressures_mpa = [0.068906986, 0.181377336, 0.424951374, 0.731306858, 1.064696014, 1.413103086], 
    concentrations_cc = [2.252715738, 5.601581157, 11.31054253, 16.87930294, 21.39238669, 25.17075548],
    temperature_k = 278.15, rho_pol_g_cm3 = 1.3937
    )
tpbo_n2_50c = IsothermData(; 
    partial_pressures_mpa = [0.269930265, 0.705742173, 1.060688385, 1.42192415, 1.813024602, 2.228349107], 
    concentrations_cc = [2.435212223, 5.677614879, 8.139676474, 10.6450967, 12.90356804, 14.82380991],
    temperature_k = 323.15, rho_pol_g_cm3 = 1.3937
    )
isotherms = [tpbo_ch4_5c, tpbo_ch4_35c, tpbo_co2_20c, tpbo_co2_50c, tpbo_n2_5c, tpbo_n2_50c]

char_co2 = [630, 300, 1.515, 44]
char_ch4 = [250, 215, 0.500, 16.04]
char_n2 = [160, 145, 0.943, 28.01]

tpbo_25_density = 1.393  # g/cm3

tpbo_co2_50c_exp = concentration(tpbo_co2_50c; component=1)
co2_bulk_phase = SL(char_co2...)

char_tpbo_valerio = [474, 900, 1.6624, 100000]
co2_tpbo25_phase_valerio = SL([474, 630], [900, 300], [1.6624, 1.515], [100000, 44], [0 0; 0 0])
co2_tpbo25_nelf_valerio = NELFModel(co2_bulk_phase, co2_tpbo25_phase_valerio, tpbo_25_density)
tpbo_co2_50c_valerio = [predict_concentration(co2_tpbo25_nelf_valerio, 323.15, p)[1] for p in partial_pressures(tpbo_co2_50c; component=1)]


#fit char params
char_tpbo25 = fit_model(NELF(), isotherms, [char_ch4, char_ch4, char_co2, char_co2, char_n2, char_n2])
@show char_tpbo25
co2_tpbo25_phase_fitted = SL([char_tpbo25[1], 630], [char_tpbo25[2], 300], [char_tpbo25[3], 1.515], [char_tpbo25[4], 44], [0 0; 0 0])
co2_tpbo25_nelf_fitted = NELFModel(co2_bulk_phase, co2_tpbo25_phase_fitted, tpbo_25_density)
tpbo_co2_50c_fitted = [predict_concentration(co2_tpbo25_nelf_fitted, 323.15, p)[1] for p in partial_pressures(tpbo_co2_50c; component=1)]

@show round.(tpbo_co2_50c_valerio; digits=2)
@show round.(tpbo_co2_50c_fitted; digits=2)
@show round.(tpbo_co2_50c_exp; digits=2)

tpbo_co2_50c_plot = plot(partial_pressures(tpbo_co2_50c; component=1), tpbo_co2_50c_valerio, label="valerio", legend=:topleft)
plot!(tpbo_co2_50c_plot, partial_pressures(tpbo_co2_50c; component=1), tpbo_co2_50c_fitted, label="fitted")
plot!(tpbo_co2_50c_plot, partial_pressures(tpbo_co2_50c; component=1), tpbo_co2_50c_exp, label="exp")