


# high res, one image


function nelf_characteristic_parameter_error_map(isotherms, bulk_phase_char_params, added_text = ""; 
    pstars=200:30:1200, tstars=200:30:1200, rhostars=missing, verbose=true, data_only=true)
    
    minimum_rhostar = maximum(polymer_density.(isotherms))
    
    if ismissing(rhostars)
        rhostars = maximum(polymer_density.(isotherms)) * 1.2
    end
    target_func = SorptionModels._make_nelf_model_parameter_target(isotherms, bulk_phase_char_params)
    needed_iters = length(pstars) * length(tstars) * length(rhostars)
    
    # start possible factored code
    errs_3D = Array{Float64, 3}(undef, length(pstars), length(tstars), length(rhostars))
    for p_i in eachindex(pstars)
        pstar = pstars[p_i]
        if verbose
            work_done = "Task " * string(round(length(tstars) * length(rhostars) * p_i / needed_iters * 100))*"% complete."
            println(work_done)
        end    
        Threads.@threads for t_i in eachindex(tstars)
            tstar = tstars[t_i]
            for r_i in eachindex(rhostars)
                rhostar = rhostars[r_i]
                errs_3D[p_i, t_i, r_i] = target_func([pstar, tstar, rhostar, DEFAULT_NELF_POLYMER_MOLECULAR_WEIGHT])
            end
        end
    end

    if typeof(rhostars) <: Number
        
        errs_2D = errs_3D[:, :, 1]
        err_2d_heatmap = contourf(
            tstars, pstars, errs_2D[:, :, 1], 
            title="rho="*string(round(rhostars[1]; digits=2))*"g/cm3" * " " * added_text, 
            xlabel = "T* (K)", ylabel = "P* (MPa)", size=(600,600), clims=(minimum(errs_2D), sum(errs_2D)/length(errs_2D)))


        if data_only; return errs_2D; end
        return err_2d_heatmap



    elseif typeof(rhostars) <: AbstractRange

        # removed factored code

        minlocation = argmin(errs_3D)
        @show global_min = [tstars[minlocation[2]]], [pstars[minlocation[1]]]
        if verbose
            println("Global minimum at a (T*, P*) of: ", global_min, " and ρ* of ", rhostars[minlocation[3]])
        end
        frameminlocations = [argmin(errs_3D[:, :, i]) for i in eachindex(rhostars)]
        mintstarvals = [tstars[frameminlocations[i][2]] for i in eachindex(rhostars)]
        minpstarvals = [pstars[frameminlocations[i][1]] for i in eachindex(rhostars)]
    

        if data_only
            return errs_3D, minlocation, frameminlocations
        end

        anim = @animate for i ∈ 1:size(errs_3D)[3]
            myplot = contourf(tstars, pstars, errs_3D[:, :, i], title="rho="*string(round(rhostars[i]; digits=3))*"g/cm3",clim=(minimum(errs_3D),5), xlabel = "T* (K)", ylabel = "P* (MPa)")
            annotate!(myplot, (0.8, 0.9), "min =" * string(round(minimum(errs_3D[:, :, i]); digits=3)))
            scatter!(myplot, [mintstarvals], [minpstarvals], legend=false)
            scatter!(myplot, global_min..., legend=false, markersize=6)
        end
        return anim
    end
   
end
