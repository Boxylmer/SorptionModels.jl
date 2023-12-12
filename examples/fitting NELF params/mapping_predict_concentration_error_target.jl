using Revise
using SorptionModels
using MembraneBase
using Plots


char_co2 = [630, 300, 1.515, 44]
char_ch4 = [250, 215, 0.500, 16.04]
char_n2 = [160, 145, 0.943, 28.01]

tpbo_25_density = 1.393  # g/cm3

co2_bulk_phase = SL(char_co2...)
# char_tpbo_valerio = [474, 900, 1.6624, 100000]
char_tpbo_valerio = [666.0811577820417, 634.7747105899552, 2.5517011086162835, 100000.0]
co2_tpbo25_phase_valerio = SL([char_tpbo_valerio[1], 630], [char_tpbo_valerio[2], 300], [char_tpbo_valerio[3], 1.515], [char_tpbo_valerio[4], 44], [0 0; 0 0])
co2_tpbo25_nelf_valerio = NELFModel(co2_bulk_phase, co2_tpbo25_phase_valerio, tpbo_25_density)




# high res, one image
pressures = 0.00:0.001:10
mass_fracs = 0.001:0.001:0.3
tpbo25_co2_error_targets = [SorptionModels._make_nelf_model_mass_fraction_target(co2_tpbo25_nelf_valerio, 308.15, p, [1]) for p in pressures]
function make_mass_frac_error_map(target_funcs, added_text = "")
    needed_iters_highres = length(pressures) * length(mass_fracs)
    errs_2D = Array{Float64, 2}(undef, length(mass_fracs), length(pressures))
    for p_i in eachindex(pressures)
        work_done = "Single image " * string(round(length(mass_fracs) * p_i / needed_iters_highres * 100))*"% complete."
        println(work_done)
        Threads.@threads for m_i in eachindex(mass_fracs)
            massfrac = mass_fracs[m_i]
            errs_2D[m_i, p_i] = target_funcs[p_i]([massfrac])
        end
    end
    # err_2d_heatmap = heatmap(tstars_highres, pstars_highres, errs_2D[:, :, 1], title="rho="*string(rhostars_highres[1])*"g/cm3", xlabel = "T* (K)", ylabel = "P* (MPa)")
    err_2d_heatmap = contourf(
        pressures, mass_fracs, errs_2D, 
        title="TPBO.25/CO2 @ 35C " * " " * added_text, 
        xlabel = "P (MPa)", ylabel = "mass frac.")
    return err_2d_heatmap
end
err_2d_heatmap = make_mass_frac_error_map(tpbo25_co2_error_targets)

expected_concs = [predict_concentration(co2_tpbo25_nelf_valerio, 308.15, p; units=:frac)[1] for p in pressures]
plot!(err_2d_heatmap, pressures, expected_concs)
# err_2d_heatmap_2 = make_mass_frac_error_map(error_target_2, "new")
# combined_plot = plot(err_2d_heatmap, err_2d_heatmap_2, layout = (1, 2))
savefig(err_2d_heatmap, joinpath(@__DIR__, "TPBO25 and CO2 at 35C and 1MPa mass frac error heatmap.png")) 
