


# high res, one image


function nelf_characteristic_parameter_error_map(isotherms, bulk_phase_char_params, added_text = ""; pstars = 200:30:1200, tstars = 200:30:1200, rhostar = missing)
    if ismissing(rhostar)
        rhostar = maximum(polymer_density.(isotherms)) * 1.2
    end
    target_func = SorptionModels._make_nelf_model_parameter_target(isotherms, bulk_phase_char_params)
    needed_iters_highres = length(pstars) * length(tstars) * length(rhostar)
    errs_2D = Array{Float64, 2}(undef, length(pstars), length(tstars))
    for p_i in eachindex(pstars)
        pstar = pstars[p_i]
        work_done = "Single image " * string(round(length(tstars) * length(rhostar) * p_i / needed_iters_highres * 100))*"% complete."
        println(work_done)
        # Threads.@threads 
        Threads.@threads for t_i in eachindex(tstars)
            tstar = tstars[t_i]
            rhostar = rhostar[1]
            errs_2D[p_i, t_i] = target_func([pstar, tstar, rhostar, DEFAULT_NELF_POLYMER_MOLECULAR_WEIGHT])
        end
    end
    # err_2d_heatmap = heatmap(tstars_highres, pstars_highres, errs_2D[:, :, 1], title="rho="*string(rhostars_highres[1])*"g/cm3", xlabel = "T* (K)", ylabel = "P* (MPa)")
    @show errs_2D
    err_2d_heatmap = contourf(
        tstars, pstars, errs_2D[:, :, 1], 
        title="rho="*string(round(rhostar[1]; digits=2))*"g/cm3" * " " * added_text, 
        xlabel = "T* (K)", ylabel = "P* (MPa)", size=(600,600), clims=(minimum(errs_2D), sum(errs_2D)/length(errs_2D)))
    return err_2d_heatmap
end



#test

