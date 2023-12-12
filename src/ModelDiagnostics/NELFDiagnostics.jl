
function nelf_characteristic_parameter_error_map(isotherms, bulk_phase_char_params, added_text = ""; 
    pstars=200:30:1200, tstars=200:30:1200, rhostars=missing, verbose=true, data_only=true, adjust_kij=false)
    if ismissing(rhostars)
        rhostars = maximum(polymer_density.(isotherms)) * 1.2
    end
    target_func = SorptionModels._make_nelf_model_parameter_target(isotherms, bulk_phase_char_params; adjust_kij)
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
    minlocation = argmin(errs_3D)
    global_min = [tstars[minlocation[2]]], [pstars[minlocation[1]]]
    if verbose
        println("Global minimum at a (T*, P*) of: ", global_min, " and ρ* of ", rhostars[minlocation[3]])
    end

    if typeof(rhostars) <: Number
        errs_2D = errs_3D[:, :, 1]
        return errs_2D
    
    elseif typeof(rhostars) <: AbstractRange
        return errs_3D
        
    end
   
end
