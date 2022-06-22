struct ZimmLundbergAnalysis{AT, AOPDT, AOPDT, CFT, ACST}
    activities::AT
    vol_fracs::VFT
    a_over_phi_values::AOPDT
    a_over_phi_derivatives::AOPDT
    cluster_functions::CFT
    average_cluster_size::ACST
end

function ZimmLundbergAnalysis(sorption_model::GABModel, activities, pen_molar_volume)
    # assumptions
    function get_volfrac(activity)
        conc = predict_concentration(sorption_model, activity)  # cc/cc
        pen_volume = conc * pen_molar_volume / MembraneBase.CC_PER_MOL_STP  # ccstp/ccpol * cm3/mol / (ccstp/mol) = (cm3/mol * mol)/ccpol = cm3/ccpol = unitless
        vol_frac = pen_volume / (1 + pen_volume)  # unitless
        return vol_frac
    end
    
    function act_per_volfrac(activity)
        vol_frac = get_volfrac(activity)
        return activity/vol_frac
    end
    
    vol_fracs = get_volfrac.(activities)
    aoverφ_vals = act_per_volfrac.(activities)
    d_aoverφ_d_a_vals = ForwardDiff.derivative.(act_per_volfrac, activities)

    # the clustering function is G_11 / molar volume
    cluster_functions = (vol_fracs .- 1) .* d_aoverφ_d_a_vals .- 1
    average_cluster_sizes = cluster_functions .* vol_fracs .+ 1

    analysis = ZimmLundbergAnalysis(activities, vol_fracs, aoverφ_vals, d_aoverφ_d_a_vals, cluster_functions, average_cluster_sizes)
    return analysis
end