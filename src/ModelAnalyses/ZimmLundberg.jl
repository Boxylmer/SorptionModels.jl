struct ZimmLundbergAnalysis{}
    average_cluster_size

end

function ZimmLundbergAnalysis(sorption_model::GABModel, activities, pen_molar_volume)
    # assumptions
    concs = sorption_model.(activities)
    pen_volumes = concs .* pen_molar_volume ./ CC_PER_MOL_STP  # ccstp/ccpol * cm3/mol / (ccstp/mol) = (cm3/mol * mol)/ccpol = cm3/pol
    vol_fracs = pen_volumes ./ (1 .+ pen_volumes)

    function act_per_volfrac(activity)
        conc = sorption_model(activity)
        pen_volume = conc * pen_molar_volume / CC_PER_MOL_STP
        vol_frac = pen_volume / (1 + pen_volume)
        return activity/vol_frac
    end

    @show d_aoverφ_d_a_vals = ForwardDiff.derivative.(act_per_volfrac, activities)

    # the clustering function is G_11 / molar volume
    @show cluster_functions = (vol_fracs .- 1) .* d_aoverφ_d_a_vals .- 1
    @show average_cluster_sizes = cluster_functions .* vol_fracs .+ 1
    return
end