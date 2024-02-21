using SorptionModels
using Test
using Measurements
using MembraneBase
using MembraneEOS
using Plots
using BenchmarkTools
using Clapeyron


precision = 5

@testset "SorptionModels.jl" begin
    # define some working isotherms
    include("tpbo_25_isotherms.jl")
    
    # test empirical models
    @testset "Empirical Sorption Models" begin
    # dual mode
        
        # single component isotherm test
        isofit = IsothermData(; partial_pressures_mpa = [0, 0.03, 0.15, 0.6, 0.9, 1.2, 1.5], concentrations_cc = [0, 1, 3, 8, 10, 12.2, 14], fugacities_mpa = [0, 0.01, 0.1, 0.6, 0.5, 1.1, 1.3])
        dmfittings = fit_dualmode_model(isofit)
        dmfittings = fit_model(DualMode(), isofit)
        dmfittings_fugacity = fit_dualmode_model(isofit, use_fugacity=true)
        
        @test dmfittings_fugacity.b != dmfittings.b

        @test round(dmfittings.ch; digits=precision) == round(6.807641216606124; digits=precision) &&
            round(dmfittings.b; digits=precision) == round(3.376206925261848; digits=precision) && 
            round(dmfittings.kd; digits=precision) == round(5.556359995999332; digits=precision)


        # assess ability to handle measurements
        isofit = IsothermData(
            partial_pressures_mpa = [0 ± 0.1, 0.03 ± 0.2, 0.15 ± 0.1, 0.6 ± 0.1, 0.9 ± 0.1, 1.2 ± 0.1, 1.5 ± 0.1], 
            concentrations_cc = [0 ± 0.3, 1 ± 0.3, 3 ± 0.3, 8 ± 0.3, 10 ± 0.3, 12.2 ± 0.3, 14 ± 0.3]
        )
        isofit_fugacity = IsothermData(
            fugacities_mpa = [0 ± 0.1, 0.03 ± 0.2, 0.15 ± 0.1, 0.6 ± 0.1, 0.9 ± 0.1, 1.2 ± 0.1, 1.5 ± 0.1], 
            concentrations_cc = [0 ± 0.3, 1 ± 0.3, 3 ± 0.3, 8 ± 0.3, 10 ± 0.3, 12.2 ± 0.3, 14 ± 0.3]
        )
        dmfittings = fit_dualmode_model(isofit)
        
        @test round(dmfittings.ch; digits=precision) == round(6.807641216606124; digits=precision) &&
            round(dmfittings.b; digits=precision) == round(3.376206925261848; digits=precision) && 
            round(dmfittings.kd; digits=precision) == round(5.556359995999332; digits=precision)


        # test jackknife method (same values as the first dmfittings)
        dmfittings_3 = fit_dualmode_model(isofit; uncertainty_method=:JackKnife)
        dmfittings_4 = fit_dualmode_model(isofit_fugacity; uncertainty_method=:JackKnife, use_fugacity=true)
        
        # print(dmfittings_3.ch.err, dmfittings_3.b.err, dmfittings_3.kd.err) (when you need to see the actual values for debugging)

        @test round(dmfittings_3.ch.err; digits=precision) == round(0.7561917289285796; digits=precision) &&
        round(dmfittings_3.b.err; digits=precision) == round(1.036574283136297; digits=precision) && 
        round(dmfittings_3.kd.err; digits=precision) == round(0.3217272537141366; digits=precision)
        # println(dmfittings_3.ch.err, ", ", dmfittings_3.b.err, ", ", dmfittings_3.kd.err) # (when you need to see the actual values for debugging)
        
        # test bootstrap method (same values as the first dmfitting) (incomplete, untested)

        # dmfittings_3 = fit_dualmode_model(isofit; uncertainty_method=:Bootstrap)
        # @test dmfittings_3.ch.err == 0.756291283433861 &&
        # dmfittings_3.b.err == 1.0368748268249859 && 
        # dmfittings_3.kd.err == 0.3217790802072902


        model = SorptionModels.DualModeModel(2.3, 4.5, 0.1)
        pressure = 5.2
        concentration = predict_concentration(model, pressure)
        recalculated_pressure = SorptionModels.predict_pressure(model, concentration)
        @test pressure == round(recalculated_pressure; digits=precision)

        # test cases
        # pentane in poly(PFMMD-co-CTFE) at 25C
        partial_pressures = [0.000536313 ± 1.34078E-06, 0.005384973 ± 1.34624E-05, 0.011818656 ± 2.95466E-05, 0.01855287 ± 4.63822E-05, 0.025361536 ± 6.34038E-05, 0.03215868 ± 8.03967E-05]
        concentrations = [0.025431809 ± 0.013099647, 1.080866801 ± 0.14091057, 2.167767053 ± 0.316078227, 3.050610392 ± 0.521260816, 3.942765687 ± 0.754846261, 4.861759146 ± 1.013138368]
        iso_case_1 = IsothermData(
            partial_pressures_mpa = partial_pressures,
            concentrations_cc = concentrations
        )
        model = fit_dualmode_model(iso_case_1)
        model_same = fit_dualmode_model(strip_measurement_to_value(iso_case_1))
        @test model.ch ≈ 1.7078301291320859
        @test model.ch == model_same.ch

        model_different = fit_dualmode_model(iso_case_1; apply_weights=true)
        # @test model.ch != model_different.ch  # todo

        # very low-solubility isotherms tend not to fit for some reason
        partial_pressures = [0.00011198, 0.000247609, 0.000521334, 0.000821977, 0.00114812]
        concentrations = [0.023724367, 0.040262462, 0.060287035, 0.072594922, 0.079677023]
        iso_low_conc = IsothermData(partial_pressures_mpa=partial_pressures, concentrations_cc = concentrations)
        model_low_conc = fit_model(DualMode(), iso_low_conc)
        @test round(model_low_conc.ch; digits=5) ≈ 0.10277
        @test round(model_low_conc.b; digits=1) ≈ 2582.4
        @test round(model_low_conc.kd; digits=3) ≈ 2.696#24


    # GAB
        acts =  [0, 0.023, 0.083, 0.161, 0.202,  0.263, 0.327 ]  # this is real 1-propanol sorption in TPBO-0.25 at 25C
        concs = [0, 2.253, 3.925, 6.126, 8.091, 12.363, 20.261]


        # test basic fitting consistency 
        gmfitting = fit_gab_model(acts, concs)
        @test gmfitting.cp ≈ 4.751529121566907
        @test gmfitting.k ≈ 2.3653742729960796
        @test gmfitting.a ≈ 9.807243268027493
        
        init_params = SorptionModels._get_initial_conditioned_gab_params(acts, concs)
        @test init_params[1] == 20.261
        @test init_params[2] ≈ 1.5758399300659984
        @test init_params[3] ≈ 0.8545460124098218



        # test checking for DimensionMismatch
        @test_throws DimensionMismatch fit_gab_model([1, 2, 3], [1, 2, 3, 4])

        # test fitting with Measurement types
        acts_2 = [0 ± 0.1, 0.023 ± 0.1, 0.083 ± 0.1, 0.161 ± 0.1, 0.202 ± 0.1, 0.263 ± 0.1, 0.327 ± 0.1]  # this is real 1-propanol sorption in TPBO-0.25 at 25C
        concs_2 = [0 ± 0.5, 2.253 ± 0.5, 3.925 ± 0.5, 6.126 ± 0.5, 8.091 ± 0.5, 12.363 ± 0.5, 20.261 ± 0.5]
        
        gmfitting_2 = fit_gab_model(acts_2, concs_2)
        @test round(gmfitting_2.cp; digits = 5) == round(4.751529121566907; digits = 5)
        @test round(gmfitting_2.k; digits = 5) ==  round(2.3653742729960796; digits = 5)
        @test round(gmfitting_2.a; digits = 5) == round(9.807243268027493; digits = 5)
        

        # test fitting with JackKnife uncertainty_method
        gmfitting_3 = fit_gab_model(acts_2, concs_2; uncertainty_method=:JackKnife)
        @test gmfitting_3.cp.err ≈ 0.4973664635441895
        @test gmfitting_3.k.err ≈ 0.12367760003383674
        @test gmfitting_3.a.err ≈ 4.090844315192642

        # very low-solubility isotherms tend not to fit for some reason
        partial_pressures = [0.00011198, 0.000247609, 0.000521334, 0.000821977, 0.00114812]
        acts = [0.035361112, 0.078190268, 0.164627403, 0.25956501, 0.362554938]
        concentrations = [0.023724367, 0.040262462, 0.060287035, 0.072594922, 0.079677023]
        iso_low_conc = IsothermData(partial_pressures_mpa=partial_pressures, concentrations_cc=concentrations, activities=acts)
        model_low_conc = fit_model(GAB(), iso_low_conc)
        init_params_low_conc = SorptionModels._get_initial_conditioned_gab_params(acts, concentrations)
        @test round(model_low_conc.cp; digits = 5) == round(0.10571589892661197; digits = 5)
        @test round(model_low_conc.k; digits = 5) == round(0.038798800644367075; digits = 5)
        @test round(model_low_conc.a; digits = 1) == 204.7#0775
        @test init_params_low_conc[1] == concentrations[end]
        @test round(init_params_low_conc[2]; digits = 5) == round(0.47376105315526396; digits = 5)
        @test round(init_params_low_conc[3]; digits = 5) == round(25.561581064859528; digits = 5)      
        # test predict_concentration using the data above 
        # easy, ideal conversion between pressure and activity
        pvap_est = sum(partial_pressures ./ acts) / length(acts)
        pressure_conversion_function(pressure) = pressure / pvap_est
        activity_conversion_function(activity) = activity * pvap_est
        gab_model_with_converters = fit_model(GAB(), iso_low_conc; pressure_conversion_function, activity_conversion_function)
        gab_model_with_converters_check = GABModel(model_low_conc.cp, model_low_conc.k, model_low_conc.a, pressure_conversion_function, activity_conversion_function)
        @test gab_model_with_converters.cp == gab_model_with_converters_check.cp
        @test round(predict_concentration(gab_model_with_converters, partial_pressures[3]); digits=5) ==
              round(a_predict_concentration(gab_model_with_converters, acts[3]); digits=5)

        # test predict activity
        given_act = acts[1]
        given_conc = a_predict_concentration(gab_model_with_converters, acts[1])
        recovered_act = predict_activity(gab_model_with_converters, given_conc)
        @test given_act ≈ recovered_act
        given_pres = partial_pressures[1]
        recovered_pres = predict_pressure(gab_model_with_converters, given_conc)
        @test round(recovered_pres; digits=8) == round(given_pres; digits=8)


        # # methanol 25C was a weird case that fit with a naieve solution but not the "sophisticated" one. This should always work. 
        # acts = []
        # concs = []
        # meth_iso = IsothermData(; activities=acts, concentrations_cc=concs)
        # gabmodel = fit_model(GAB(), meth_iso)

    
        # Flory Huggins - Dual Mode
    #   TPBO-1.00 with propanol at 25C
    acts_3 = [0.01061664, 0.034878645, 0.06508188, 0.091689136, 0.126116896, 0.163921419]
    concs_3 = [1.188187617 ± 0.013431227, 2.455888283 ± 0.027591385, 3.80430418 ± 0.042799107, 5.824894435 ± 0.06532282, 9.604100032 ± 0.107154714, 14.70392798 ± 0.163320976]
    iso_3 = IsothermData(; activities=acts_3, concentrations_cc = concs_3)
    fhdm = fit_model(FloryHugginsDualMode(), iso_3, 74.7453)
    pred = [a_predict_concentration(fhdm, a) for a in acts_3]

    end

    @testset "Fundamental Sorption Models" begin
        using MembraneEOS
        polymer = "PC"
        penetrant = "CO2"
        kij = [0 -0.007; -0.007 0]
        bulk_phase_eos = MembraneEOS.SL([penetrant])
        polymer_phase_eos = MembraneEOS.SL([polymer, penetrant], kij)
        density = 1.197850471   #g/cm3    
        temperature = 308.15
        pressures = [0, 0.18, 0.38, 0.64, 0.94, 1.23, 1.44]
        expected_mass_fracs = [5.51E-06, 0.008294923, 0.014447025, 0.020467468, 0.026066798, 0.030734723, 0.033827052]
        expected_concs_cc_cc = [0, 5.094493596, 8.910071012, 12.66689539, 16.17497812, 19.10613428, 21.05001223]
        
        # NELF    
            ksw = 0.0102            # 1/MPa
            nelfmodel = NELFModel(bulk_phase_eos, polymer_phase_eos, density)
            nelf_concs_pure_co2 = [predict_concentration(nelfmodel, temperature, p, [1.0]; ksw=[ksw])[1] for p in pressures]
        
        
            penetrants = ["CO2", "CO2"]
            kij_ternary = [0      -0.007 -0.007; 
                           -0.007 0      0.0   ; 
                           -0.007 0.0    0.0   ]
            ksw_ternary = [0.0102, 0.0]            # 1/MPa
            bulk_phase_eos_ternary = MembraneEOS.SL(penetrants)
            polymer_phase_eos_ternary = MembraneEOS.SL([polymer, penetrants...], kij_ternary)
            nelfmodel_ternary = NELFModel(bulk_phase_eos_ternary, polymer_phase_eos_ternary, density)
            nelf_concs_co2_mix = [predict_concentration(nelfmodel_ternary, temperature, p, [0.5, 0.5]; ksw=ksw_ternary)[1] for p in pressures]
            @test nelf_concs_co2_mix[3] != nelf_concs_pure_co2[3]
        
            nelf_concs_co2_psuedo = [predict_concentration(nelfmodel_ternary, temperature, p, [1.0, 0]; ksw=ksw_ternary)[1] for p in pressures]
            @test nelf_concs_co2_psuedo[3] ≈ nelf_concs_pure_co2[3]

            @test round(infinite_dilution_solubility(nelfmodel, temperature)) ≈ 40

            # test the polymer fitter with TPBO-0.25
            char_co2 = [630, 300, 1.515, 44]
            char_ch4 = [250, 215, 0.500, 16.04]
            char_n2 = [160, 145, 0.943, 28.01]

            # char_tpbo_valerio = [474, 900, 1.6624, 100000]

            isotherms = [tpbo_ch4_5c, tpbo_ch4_20c, tpbo_ch4_35c, tpbo_co2_5c, tpbo_co2_20c, tpbo_co2_35c, tpbo_co2_50c, tpbo_n2_5c, tpbo_n2_50c]
            bulk_phase_char_params = [char_ch4, char_ch4, char_ch4, char_co2, char_co2, char_co2, char_co2, char_n2, char_n2]
            char_tpbo25 = fit_model(NELF(), isotherms, bulk_phase_char_params, verbose=false; initial_search_resolution=10)

            # just running to make sure it doesn't throw
            char_tpbo25_with_errs_hessian = fit_model(NELF(), isotherms, bulk_phase_char_params, verbose=false; initial_search_resolution=10, uncertainty_method=:Hessian)
            
            
            # infinite dilution
            sinf_entr = infinite_dilution_solubility_entropic(nelfmodel)
            sinf_enth = infinite_dilution_solubility_enthalpic(nelfmodel, temperature)

            ln_sinf = sinf_enth + sinf_entr + log(273.15/(0.101325 * temperature))
            sinf_analytical = exp(ln_sinf)
            sinf_numerical = infinite_dilution_solubility(nelfmodel, temperature)
            @test round(sinf_analytical; sigdigits=2) == round(sinf_numerical; sigdigits=2)
            
            
            
            
            # dualmode_models = [fit_model(DualMode(), isotherm) for isotherm in isotherms]
            # now that we have some characteristic parameters, we can try to fit individual kij and ksw for a gas
            #   pick CO2
            # given_co2_50c = concentration(tpbo_co2_50c; component=1)
            # given_co2_20c = concentration(tpbo_co2_20c; component=1)
            
            # co2_bulk_phase = SL(char_co2...)
            
            # polymer_phase_valerio = SL([474, 630], [900, 300], [1.6624, 1.515], [100000, 44], [0 -0.0356; -0.0356 0])
            # nelf_model_valerio = NELFModel(co2_bulk_phase, polymer_phase_valerio, 1.393)
            # given_valerio_co2_50c = [predict_concentration(nelf_model_valerio, 323.15, p)[1] for p in partial_pressures(tpbo_co2_50c; component=1)]
            # given_valerio_co2_20c = [predict_concentration(nelf_model_valerio, 293.15, p)[1] for p in partial_pressures(tpbo_co2_20c; component=1)]
            # @show fit_kij(NELF(), [tpbo_co2_20c, tpbo_co2_50c], char_co2, char_tpbo_valerio)
            # @show fit_ksw(NELF(), tpbo_co2_20c, co2_bulk_phase, polymer_phase_valerio)
            # fit char params to exp data
            # char_tpbo25 = fit_model(NELF(), isotherms, [char_ch4, char_ch4, char_co2, char_co2, char_n2, char_n2])
            # kij_co2_tpbo25 = fit_kij(NELF(), [tpbo_co2_20c, tpbo_co2_50c], char_co2, char_tpbo25)
            # polymer_phase_fit_no_kij = SL([char_tpbo25[1], 630], [char_tpbo25[2], 300], [char_tpbo25[3], 1.515], [char_tpbo25[4], 44], [0 0; 0 0])
            # polymer_phase_fit_with_kij = SL([char_tpbo25[1], 630], [char_tpbo25[2], 300], [char_tpbo25[3], 1.515], [char_tpbo25[4], 44], [0 kij_co2_tpbo25; kij_co2_tpbo25 0])
            # ksw_co2_tpbo_20c = fit_ksw(NELF(), tpbo_co2_20c, co2_bulk_phase, polymer_phase_fit_with_kij)
            # fit_pred_with_kij_and_ksw_co2_20c = [predict_concentration(nelf_model_fit_with_kij, 293.15, p; ksw=ksw_co2_tpbo_20c)[1] for p in partial_pressures(tpbo_co2_20c; component=1)]
            
            
            
            
            # make predictions
            # nelf_model_fit = NELFModel(co2_bulk_phase, polymer_phase_fit_no_kij, 1.393)
            # nelf_model_fit_with_kij = NELFModel(co2_bulk_phase, polymer_phase_fit_with_kij, 1.393)
            # nelf_model_fit_with_kij = NELFModel(co2_bulk_phase, polymer_phase_fit_with_kij, 1.393)

            # fit_pred_no_kij_co2_50c = [predict_concentration(nelf_model_fit, 323.15, p)[1] for p in partial_pressures(tpbo_co2_50c; component=1)]
            # fit_pred_with_kij_co2_50c = [predict_concentration(nelf_model_fit_with_kij, 323.15, p, [1])[1] for p in partial_pressures(tpbo_co2_50c; component=1)]
            
            # fit_pred_no_kij_co2_20c = [predict_concentration(nelf_model_fit, 293.15, p)[1] for p in partial_pressures(tpbo_co2_20c; component=1)]
            # fit_pred_with_kij_co2_20c = [predict_concentration(nelf_model_fit_with_kij, 293.15, p)[1] for p in partial_pressures(tpbo_co2_20c; component=1)]


            # @show char_tpbo25
            # @show kij_co2_tpbo25
            # @show ksw_co2_tpbo_20c
            # @show round.(given_co2_50c)
            # @show round.(given_valerio_co2_50c)
            # @show round.(fit_pred_no_kij_co2_50c)
            # @show round.(fit_pred_with_kij_co2_50c)
            # @show ""
            # @show round.(given_co2_20c)
            # @show round.(given_valerio_co2_20c)
            # @show round.(fit_pred_no_kij_co2_20c)
            # @show round.(fit_pred_with_kij_co2_20c)
            # @show round.(fit_pred_with_kij_and_ksw_co2_20c)


            # co2_sl = SL("CO2")
            # pc_sl = SL(["PC", "CO2"])
            # co2_pc_nelf = NELFModel(co2_sl, pc_sl, 1.197850471)
            # predict_concentration(co2_pc_nelf, 308.15, 0.38, [1]; ksw = [0.0102])[1]   # should be 8.91 to 8.92
            # ps = [1.00E-05, 0.157894737, 0.315789474, 0.473684211, 0.631578947, 0.789473684, 0.947368421, 1.105263158, 1.263157895, 1.421052632, 1.578947368, 1.736842105, 1.894736842, 2.052631579, 2.210526316, 2.368421053, 2.526315789, 2.684210526, 2.842105263, 3, ]
            # expected_res = [0.003361992, 4.365419237, 7.136022736, 9.144646417, 10.70671881, 11.97602153, 13.03904388, 13.94920076, 14.74174998, 15.44115597, 16.06506788, 16.62662242, 17.13585174, 17.60058337, 18.02703695, 18.42023244, 18.78427702, 19.12257132, 19.43796081, 19.73284857, ]
            # res = [predict_concentration(co2_pc_nelf, 308.15, p, [1]; ksw = [0.0])[1] for p in ps]  # no swelling or kij
            
        
        # DGRPT

            dgrptmodel = DGRPTModel(bulk_phase_eos, polymer_phase_eos, density)
            # dgrpt_concs_pure_co2 = [predict_concentration(dgrptmodel, temperature, p, [1.0])[1] for p in pressures]

            # compare to nelf w/ no swelling
            # ps = [0.1, 1, 2, 3, 3.5]
            ps = [0.1, 3.5]
            nelfmodel = NELFModel(bulk_phase_eos, polymer_phase_eos, density)
            
            # nelf_concs_pure_co2 = [predict_concentration(nelfmodel, temperature, p, [1.0])[1] for p in ps]
            # nelf_densities = [calculate_swelled_polymer_density(nelfmodel, p, [1], [0]) for p in ps]
            # dgrpt_concs_pure_co2_1_term = [predict_concentration(dgrptmodel, temperature, p, [1.0]; taylor_series_order=1)[1] for p in ps]
            # dgrpt_dens_pure_co2_1_term = [polymer_density(dgrptmodel, temperature, p, [1]; taylor_series_order=1) for p in ps]
            # dgrpt_concs_pure_co2_2_terms = [predict_concentration(dgrptmodel, temperature, p, [1.0]; taylor_series_order=2)[1] for p in ps]
            # dgrpt_concs_pure_co2_5_terms = [predict_concentration(dgrptmodel, temperature, p, [1.0]; taylor_series_order=5)[1] for p in ps]
            # @show nelf_concs_pure_co2
            # @show nelf_densities
            # @show dgrpt_concs_pure_co2_1_term
            # @show dgrpt_dens_pure_co2_1_term
            # @show dgrpt_concs_pure_co2_2_terms
            # @show dgrpt_concs_pure_co2_5_terms
    end

    @testset "Transient Sorption Models" begin
        
        fickian_model = FickianSorptionModel(1, 0.01)
        bh_model = BerensHopfenbergSorptionModel(0.7, 0.01, 0.3, 0.001)
        mbh_model = ModifiedBerensHopfenbergSorptionModel(0.7, 0.01, 0.3, 0.001, 0.05)
        fickian_prediction = predict_sorption(fickian_model, 10)
        bh_prediction = predict_sorption(bh_model, 10)
        mbh_prediction = predict_sorption(mbh_model, 10)

        time_vector = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30, 40, 50, 70, 90, 100]
        fick_prediction_vector = predict_sorption(fickian_model, time_vector)
        bh_prediction_vector = predict_sorption(bh_model, time_vector)
        mbh_prediction_vector = predict_sorption(mbh_model, time_vector)
        
        
        # testing fitting without errors (no true tests)
        time_data   =        [1.,       2,       3,        4,       5,      6,      7,      8,      9,      10,    15,    20,    25,    30,    40,   50,    70,   90,     100] 
        transient_sorption = [0.00194, 0.00515, 0.009107, 0.01359, 0.0184, 0.0237, 0.0291, 0.0348, 0.0407, 0.046, 0.077, 0.109, 0.141, 0.171, 0.22, 0.278, 0.36, 0.4364, 0.4671]   
        exp_data = TransientStepData(time_data, transient_sorption)

        fick_model_fit = fit_transient_sorption_model(exp_data, FickianSorption())
        
        bh_model_fit = fit_transient_sorption_model(exp_data, BerensHopfenbergSorption())
    
        mbh_model_fit = fit_transient_sorption_model(exp_data, ModifiedBerensHopfenbergSorption())
        
        semi_thickness_cm = 0.02 ± 0.001

        fick_d = get_diffusivity(fick_model_fit, semi_thickness_cm)
        @test isapprox(fick_d.val, 4.779020339107122e-7, atol =1e-5)
        @test isapprox(fick_d.err, 4.7790203391071216e-8, atol=1e-7)

        mbh_d = get_diffusivity(mbh_model_fit, semi_thickness_cm)
        @test isapprox(mbh_d.val, 2.251491971659755e-15, atol=1e-14)
        @test isapprox(mbh_d.err, 2.251491971659755e-16, atol=1e-15)

        # boostrap and uncertainty
        fick_model_fit_with_err = fit_transient_sorption_model(exp_data, FickianSorption(); uncertainty_method=:Bootstrap)
        fick_d_with_err = get_diffusivity(fick_model_fit_with_err, semi_thickness_cm)
        @test fick_d_with_err.val ≈ fick_d.val
        # @test fick_d_with_err.err ≈ 7.950417186869458e-8 # errors are random due to low sample size
        
        shifted_exp_data = TransientStepData(time_data .+ 20., transient_sorption)
        sbh_model_fit = fit_transient_sorption_model(shifted_exp_data, ModifiedBerensHopfenbergSorption())
        sbh_prediction_vector = predict_sorption(sbh_model_fit, shifted_exp_data.time)


        # benchmarks and allocating
        #10,546,944
        #11,096,576
        #272,784
        single_fit_allocs = @allocated fit_transient_sorption_model(exp_data, FickianSorption())
        #129,185,280
        #140,445,520
        #132,250,512
        #3,434,224
        #3,341,584
        bootstrap_fit_allocs = @allocated  fit_transient_sorption_model(exp_data, FickianSorption(); uncertainty_method=:Bootstrap)
        # @test bootstrap_fit_allocs <= 3341584
        #142.868 ms
        #19.243 ms
        # @show @btime fit_transient_sorption_model($exp_data, BerensHopfenbergSorption(); uncertainty_method=:Bootstrap)
        # @profview fit_transient_sorption_model(exp_data, FickianSorption(); uncertainty_method=:Bootstrap)
    end

    @testset "Dilation Models" begin
        pressures_mpa = ([0.0, 4.34, 9.04, 12.8, 17.0, 21.16, 25.15] .± 0.01) .* 0.101325 
        frac_dilations = ([0.0, 0.4778, 1.01, 1.458, 1.93, 2.44, 2.93] .± 0.02) ./ 100 
        
        # empirical dilation model
        edm = fit_model(EmpiricalDilation(), pressures_mpa, frac_dilations)
        edm_jk = fit_model(EmpiricalDilation(), pressures_mpa, frac_dilations, :JackKnife)
        edm_h = fit_model(EmpiricalDilation(), pressures_mpa, frac_dilations, :Hessian)
        
        edm_n1 = fit_model(EmpiricalDilation(), pressures_mpa, frac_dilations, n_params = 1)
        edm_n2 = fit_model(EmpiricalDilation(), pressures_mpa, frac_dilations, n_params = 2)
        edm_n3 = fit_model(EmpiricalDilation(), pressures_mpa, frac_dilations, n_params = 3)
        edm_n4 = fit_model(EmpiricalDilation(), pressures_mpa, frac_dilations, n_params = 4)

        pred_edm = predict_dilation(edm, pressures_mpa)
        pred_edm_jk = predict_dilation(edm_jk, pressures_mpa)
        pred_edm_h = predict_dilation(edm_h, pressures_mpa)

        @test edm == edm_n3

        @test strip_measurement_to_value(pred_edm) == strip_measurement_to_value(pred_edm_jk) == strip_measurement_to_value(pred_edm_h)
        @test edm_n1 != edm_n2 != edm_n3 != edm_n4
       
        @test size(predict_dilation_derivative(edm, pressures_mpa))[1] == 7
        
        # dual mode dilation model
        dmmodel = DualModeModel(0.01, 0.01, 7.8037/(16.698 * 0.101325) ± 0.05) # cc/mpa
        dmmodel_throws = DualModeModel(0.01, 0.01, 7.8037/(16.698 * 0.101325) ± 0.05; use_fugacity=true) # cc/mpa
        @test_throws ArgumentError fit_model(DualModeDilation(), pressures_mpa, frac_dilations, dmmodel_throws)
        dmd = fit_model(DualModeDilation(), pressures_mpa, frac_dilations, dmmodel)
        dmd_jk = fit_model(DualModeDilation(), pressures_mpa, frac_dilations, dmmodel, :JackKnife)
        dmd_h = fit_model(DualModeDilation(), pressures_mpa, frac_dilations, dmmodel, :Hessian)
        @test_throws ArgumentError fit_model(DualModeDilation(), pressures_mpa, frac_dilations, dmmodel, :InvalidErrorProp)
        pred_dmd = predict_dilation(dmd, pressures_mpa)
        pred_dmd_jk = predict_dilation(dmd_jk, pressures_mpa)
        pred_dmd_h  = predict_dilation(dmd_h, pressures_mpa)

        @test strip_measurement_to_value(pred_dmd) == strip_measurement_to_value(pred_dmd_jk) == strip_measurement_to_value(pred_dmd_h)
    
        struct BadDilationModel <: SorptionModels.DilationModel end
        unimplemented_model = BadDilationModel()
        @test_throws ErrorException predict_dilation(unimplemented_model)
        # @test_throws f
        
    end

    isotherm_1 = IsothermData(;  # CH4 in TPBO-0.50 at 20C
        partial_pressures_mpa = [0.259163083, 0.680732225, 1.153210606, 1.681231582, 2.230997679, 2.726364496, 3.143840263],
        concentrations_cc = [20.06902712, 34.42914464, 44.07950458, 51.56463956, 57.01452812, 61.02024241, 63.74277795],
        temperature_k = 293.15,
        fugacities_mpa = [0.257596406, 0.670008659, 1.122709817, 1.617057535, 2.119186153, 2.560995079, 2.925750142]
    )
    isotherm_2 = IsothermData(; 
        partial_pressures_mpa = [0.262492618, 0.683659348, 1.155742537, 1.686317619, 2.20405451, 2.706619691, 3.136376749],
        concentrations_cc = [17.57118063, 31.09445449, 39.7469948, 46.41292601, 51.16825653, 54.8539774, 56.54801589],
        temperature_k = 300.15,
        fugacities_mpa = [0.261005294, 0.67364989, 1.127391615, 1.626570478, 2.103003269, 2.555717054, 2.935454538]
    )
    isotherm_3 = IsothermData(; 
        partial_pressures_mpa = [0.266511453, 0.686673842, 1.157709902, 1.687453932, 2.2030762, 2.842052813, 3.197348035],
        concentrations_cc = [15.05920737, 26.87379577, 34.86510363, 40.73884457, 45.15400583, 49.76327481, 51.51256356],
        temperature_k = 308.15,
        fugacities_mpa =[0.265107381, 0.677426321, 1.131657671, 1.632663035, 2.110610706, 2.690078663, 3.006342365]
    )
    isotherm_4 = IsothermData(; 
        partial_pressures_mpa = [0.274645736, 0.692967011, 1.177848036, 1.698686751, 2.207084424, 2.71751992, 3.142763943],
        concentrations_cc = [11.78645049, 22.51120434, 30.36900358, 36.63755063, 40.37998354, 44.29098419, 46.8613559],
        temperature_k = 323.15,
        fugacities_mpa =[0.27337983, 0.684971536, 1.154961289, 1.651558222, 2.128304942, 2.599274758, 2.985937811])
    isotherms = [isotherm_1, isotherm_2, isotherm_3, isotherm_4]



    g_isotherm_1 =  IsothermData(; 
        partial_pressures_mpa = [0.000138386, 0.000583942, 0.00104234, 0.001498221],
        concentrations_cc = [0.682951321, 1.111629734, 1.407396085, 1.656181954],
        activities = [0.04369976, 0.184397831, 0.329151532, 0.473110012],
        temperature_k = 298.15)
    g_isotherm_2 =  IsothermData(; 
        partial_pressures_mpa = [0.000348261, 0.000916842, 0.001463913, 0.002013374, 0.002596481, 0.003140448],
        concentrations_cc = [0.731497769, 1.1041389, 1.293588162, 1.597390536, 1.949575609, 2.422564057],
        activities = [0.061926184, 0.163028384, 0.260305964, 0.358008582, 0.461693813, 0.558419375],
        temperature_k = 308.15)
    g_isotherm_3 =  IsothermData(; 
        partial_pressures_mpa = [0.000637003, 0.001775022, 0.002649602],
        concentrations_cc = [0.919657573, 1.292284349, 1.53203915],
        activities = [0.066452474, 0.185171251, 0.276407799],
        temperature_k = 318.15)
    g_isotherm_pvaps = [0.003166749, 0.005623816, 0.009585843]
    g_isotherms = [g_isotherm_1, g_isotherm_2, g_isotherm_3]



    @testset "Analysis Functions" begin
        
        # VHDM analysis
        @testset "Van't Hoff Dual Mode Fitting" begin    
            vhdm_analysis = VantHoffDualModeAnalysis(isotherms)
            vhdm_analysis_but_with_fugacity = VantHoffDualModeAnalysis(isotherms; use_fugacity=true)
            @test vhdm_analysis_but_with_fugacity.final_models[1].use_fugacity == true
            @test vhdm_analysis_but_with_fugacity.final_models[2].ch.val ≈ 57.13799466676326
            @test round(vhdm_analysis_but_with_fugacity.final_models[2].ch.err; digits=4) ≈ 29.838700113541947
            @test vhdm_analysis.final_models[2].ch.val ≈ 56.37550904605368
        end

        # Isosteric Heat
        @testset "Isosteric Heat Analysis" begin    
            gab_pressure_conversion_funcs = [(p) -> p/pvap for pvap in g_isotherm_pvaps]
            gab_activity_conversion_funcs = [(a) -> a*pvap for pvap in g_isotherm_pvaps]
            
            eosmodel(p, t) = MembraneEOS.compressibility_factor(MembraneEOS.PR("CH4"), p, t) 
            unideal_ish_analysis = IsostericHeatAnalysis(isotherms, eosmodel) 
            ish_analysis = IsostericHeatAnalysis(isotherms)

            no_err_isos = strip_measurement_to_value.(isotherms)
            no_err_isosteric_heat = IsostericHeatAnalysis(no_err_isos, linear_regression_error=false)
            @test typeof(no_err_isosteric_heat.isosteric_heat_at_conc[1]) <: Float64
            
            @test unideal_ish_analysis.isosteric_heat_at_conc[2] ≠ ish_analysis.isosteric_heat_at_conc[2]
            ish_analysis_but_with_vhdm = IsostericHeatAnalysis(isotherms; use_vant_hoff_constraints=true)
            @test_throws ErrorException IsostericHeatAnalysis(isotherms; model=GAB())
            @test_nowarn ish_analysis_but_with_gab = IsostericHeatAnalysis(g_isotherms; model=GAB(), gab_pressure_conversion_funcs, gab_activity_conversion_funcs)
            
            conc_no_constraint = ish_analysis.sampled_concentrations[4]
            conc_with_constraint = ish_analysis_but_with_vhdm.sampled_concentrations[4]
            @test conc_no_constraint == conc_with_constraint
            @test ish_analysis.isosteric_heat_at_conc[4].val ≈ -18815.424127650764
            @test ish_analysis_but_with_vhdm.isosteric_heat_at_conc[4].val ≈ -18121.602388309566
        end

        @testset "Webb Isosteric Heat Analysis" begin
            eos_ch4 = MembraneEOS.PR("CH4")
            eos_z(p, t) = MembraneEOS.compressibility_factor(eos_ch4, p, t)
            wish_analysis_no_eos = SorptionModels.WebbIsostericHeatAnalysis(isotherms)
            wish_analysis = SorptionModels.WebbIsostericHeatAnalysis(isotherms, eos_z)

        end

        # Mobility Factor
        @testset "Mobility Factor and Thermo Factor" begin
            # water diffusion in TPBO-0.25 at 25C, taken from Box et al., in BH model 
            act = [0.067952153, 0.128795589, 0.184686884, 0.226021282, 0.261604807, 0.317356324, 0.416858405, 0.547610888, 0.640042101]
            dif = [9.54865E-07, 9.14651E-07, 8.29247E-07, 1.36793E-06, 8.94442E-07, 8.86731E-07, 8.78362E-07, 6.20711E-06, 5.29102E-06]

            con = [4.553301251, 7.360868935, 9.972461154, 12.76467679, 14.32852228, 17.19540216, 22.19415472, 28.84150735, 33.56667178]
            err = [0.04247127, 0.072640132, 0.104512929, 0.134097291, 0.161762819, 0.199449306, 0.257028442, 0.333787113, 0.403308077]
            con = [con[i] ± err[i] for i in eachindex(con, err)]

            iso = IsothermData(; 
                # partial_pressures_mpa=nothing, 
                concentrations_cc=con, 
                activities=act, 
                # temperature_k=nothing, 
                rho_pol_g_cm3=1.393 ± 0.002, 
                pen_mws_g_mol=18.01528
            )
            mob_fact_analysis = MobilityFactorAnalysis(dif, iso)
            therm_fact_analysis = ThermodynamicFactorAnalysis(iso)
            @test mob_fact_analysis.kinetic_factors[3].val ≈ 8.117024704029978e-7
            @test mob_fact_analysis.thermo_factor_analysis.thermodynamic_factors[3].val ≈ 1.0216144834304761
            @test therm_fact_analysis.thermodynamic_factors[3] == mob_fact_analysis.thermo_factor_analysis.thermodynamic_factors[3]


            # now test with actual acitivity functions (the method which does not use activity approximations)
            function penetrant_activity(
                    p_mpa::Number; 
                    T_ref = 160.15, 
                    model=Clapeyron.PR(["Methane"])
                )

                Saturation_p_ref = saturation_pressure(model, T_ref)[1]
                mu_ref = Clapeyron.chemical_potential(model, Saturation_p_ref, T_ref)[1]
                p = p_mpa * 1e6
                mu2 = Clapeyron.chemical_potential(model, p, 308.15)[1]
                a = exp((mu2-mu_ref)/(8.314*308.15))
                return a
            end

            CPIM_CH4_ISOTHERM = IsothermData(
                partial_pressures_mpa=[0.309802563, 0.710304653, 1.222477562, 1.832018444, 2.840712537],
                concentrations_cc = [8.147067101, 14.38986172, 20.52870123, 25.0722395, 32.85525989],
                rho_pol_g_cm3 = 1.285,
                pen_mws_g_mol = 16.043)

        
            CPIM_CH4_DM = fit_model(DualMode(),CPIM_CH4_ISOTHERM; uncertainty_method = :JackKnife)
            CPIM_CH4_TFA_1 = ThermodynamicFactorAnalysis(CPIM_CH4_ISOTHERM, CPIM_CH4_DM, x -> penetrant_activity(x, T_ref = 160.15))
            CPIM_CH4_TFA_2 = ThermodynamicFactorAnalysis(CPIM_CH4_ISOTHERM, CPIM_CH4_DM, x -> penetrant_activity(x, T_ref = 180.15))

            @test CPIM_CH4_TFA_1.thermodynamic_factors[2] == CPIM_CH4_TFA_2.thermodynamic_factors[2]
            @test CPIM_CH4_TFA_1.lna[2] != CPIM_CH4_TFA_2.lna[2]
            @test CPIM_CH4_TFA_1.lnw[2] == CPIM_CH4_TFA_2.lnw[2]

            
            # mobility factor analysis (direct method) 
            CPIM_CH4_Diffusivities = [6.00E-08, 7.21E-08, 8.27E-08, 9.68E-08]
            CPIM_CH4_D_Pressures = [0.108381747, 0.352283335, 0.632820204, 0.995498696]
            CPIM_CH4_L_ANALYSIS = MobilityFactorAnalysis(CPIM_CH4_Diffusivities, CPIM_CH4_D_Pressures, 1.285, 16.043, CPIM_CH4_DM, penetrant_activity)
            
            CPIM_CH4_Diffusivities = [6.00E-08 ± 1.93E-08, 7.21E-08 ± 2.07E-08, 8.27E-08 ± 2.18E-08, 9.68E-08 ± 2.37E-08] # now do it with uncertainty!
            CPIM_CH4_L_ANALYSIS = MobilityFactorAnalysis(CPIM_CH4_Diffusivities, CPIM_CH4_D_Pressures, 1.285, 16.043, CPIM_CH4_DM, penetrant_activity)

            
        end
        
        # Partial Immobilization Model
        @testset "Partial Immobilization Model" begin
            model = DualModeModel(56.8, 7.4, 26.1; use_fugacity=true)
            permeabilities = [1221, 1091, 1038]
            pressures_mpa = [0.298, 0.637, 0.974]
            result = PartialImmobilizationModel(model, pressures_mpa, permeabilities)
            @test result.henry_mode_diffusivity.val ≈ 2.649330314735546e-6
            @test result.langmuir_mode_diffusivity.val ≈ 1.7126341427211583e-7
        end

        # Zimm Lundberg Analysis
        @testset "Zimm Lundberg Analysis" begin
            # water isotherm taken from from TPBO-0.25 at 25C
            activities=[0.068087235, 0.129051621, 0.185054022, 0.226470588, 0.26212485, 0.317987195, 0.417687075, 0.54869948, 0.641314436]
            concentrations_cc = [4.553554612, 7.361278518, 9.973016055, 12.76538706, 14.32931957, 17.19635897, 22.19538967, 28.84311218, 33.56853954]
            isotherm = IsothermData(; activities, concentrations_cc)
            gabmodel = fit_model(GAB(), isotherm)
            
            mol_vol = 18 # cm3/mol
            zma1 = ZimmLundbergAnalysis(gabmodel, activities, mol_vol)

            # do errors work as well?
            mol_vol = 18 ± 0.2
            activities = activities .± (activities .* 0.001); concentrations_cc = concentrations_cc .± (concentrations_cc .* 0.001)
            isotherm = IsothermData(; activities, concentrations_cc)
            gabmodel = fit_model(GAB(), isotherm)
            zma2 = ZimmLundbergAnalysis(gabmodel, activities, mol_vol)

            # make sure nothing too terriblehappened
            @test zma1.average_cluster_size[1] == zma2.average_cluster_size[1].val 
            @test zma1.a_over_phi_derivatives[1] == zma2.a_over_phi_derivatives[1].val 
        end

        # Dual Mode Desorption
        @testset "Dual Mode Desorption Analysis" begin
            isotherm_desorption = IsothermData(; 
                partial_pressures_mpa = [0.241333352, 0.600763584, 1.04806673, 1.466095481, 1.951571285, 2.499847618, 3.142683031, 2.60974313, 1.199714642, 0.575992209, 0.30991402, 0.145032338],
                concentrations_cc = [41.48924079, 62.79313671, 77.9590348, 88.08019013, 96.61909374, 105.7659302, 114.9523482, 111.9833899, 98.93545196, 81.21343867, 66.62816092, 50.92125421]
            )
            isotherm_desorption_fugacity_bad = IsothermData(; 
                fugacities_mpa = [0.241333352, 0.600763584, 1.04806673, 1.466095481, 1.951571285, 2.499847618, 3.142683031, 2.60974313],
                concentrations_cc = [41.48924079, 62.79313671, 77.9590348, 88.08019013, 96.61909374, 105.7659302, 114.9523482, 111.9833899]
            )
            DualModeDesorption(isotherm_desorption_fugacity_bad; use_fugacity=true) # run with fugacity
            @test_logs DualModeDesorption(isotherm_desorption_fugacity_bad; use_fugacity=true, verbose=true)
            isotherm_very_bad = IsothermData(; 
                partial_pressures_mpa = [0.241333352, 0.600763584, 1.04806673, 1.466095481, 1.951571285, 2.499847618, 3.142683031],
                concentrations_cc = [41.48924079, 62.79313671, 77.9590348, 88.08019013, 96.61909374, 105.7659302, 114.9523482]
            )
            @test isnothing(DualModeDesorption(isotherm_very_bad))

            dmda = DualModeDesorption(isotherm_desorption)
            dmda_naive = DualModeDesorption(isotherm_desorption; naive=true)
            dmda_static_b = DualModeDesorption(isotherm_desorption; share_b=false)
            
            dmda_with_err = DualModeDesorption(isotherm_desorption; uncertainty_method=:JackKnife)
            dmda_naive_with_err = DualModeDesorption(isotherm_desorption; naive=true, uncertainty_method=:JackKnife)
            dmda_static_b_with_err = DualModeDesorption(isotherm_desorption; share_b=false, uncertainty_method=:JackKnife)

            @test_throws ArgumentError DualModeDesorption(isotherm_desorption; uncertainty_method=:invalid)

            @test round(dmda_with_err.sorbing_model.ch.val; digits=1) == 71.3
            @test round(dmda_with_err.sorbing_model.ch.err; digits=1) == 2.6

            @test round(dmda_static_b_with_err.desorbing_model.ch.val; digits=1) == 124.7
            @test round(dmda_static_b_with_err.desorbing_model.ch.err; digits=1) == 2.7

            # non-desorbing isotherms shouldn't work
            @test isnothing(DualModeDesorption(isotherm_1))

            # test a case where an isotherm with only one desorption step didn't work. This is because jackknife doesn't work when one missing data point causes there to be no desorption data.
            pres = [1 ± 0.1, 2, 3, 4, 3]
            conc = [5, 9 ± 0.4, 12, 17, 13]
            iso = IsothermData(partial_pressures_mpa = pres, concentrations_cc = conc)
            @test_nowarn DualModeDesorption(iso; uncertainty_method=:JackKnife)


            @test round(predict_concentration(dmda_naive_with_err, pres)[1][1].val) == 77
            @test predict_concentration(dmda_naive_with_err, iso)[1][1] == predict_concentration(dmda_naive_with_err, partial_pressures(iso, component=1))[1][1]
            # direct_pressures = partial_pressures(isotherm_desorption; component=1)
            # direct_concentrations = concentration(isotherm_desorption; component=1)
            # predicted_pressures = 0:0.05:maximum(direct_pressures)
            # predicted_sorbing, predicted_desorbing = predict_concentration(dmda, predicted_pressures)
        end

        # Partial Molar Volumes
        @testset "Molar Volume Analysis" begin
            # no uncertainty
            model = DualModeModel(0, 0, 7.8037/(16.698 * 0.101325)) # cc/mpa
            pressures_mpa = [0.0, 4.34, 9.04, 12.8, 17.0, 21.16, 25.15] .* 0.101325
            frac_dilations = [0.0, 0.4778, 1.01, 1.458, 1.93, 2.44, 2.93] ./ 100
            molar_vol_analysis = MolarVolumeAnalysis(model, pressures_mpa, frac_dilations)
            res_no_uncertain = molar_vol_analysis.partial_molar_volumes_cm3_mol[3]
    
            # with uncertainty
            model = DualModeModel(0, 0, 7.8037/(16.698 * 0.101325) ± 0.05) # cc/mpa
            pressures_mpa = ([0.0, 4.34, 9.04, 12.8, 17.0, 21.16, 25.15] .± 0.01) .* 0.101325 
            frac_dilations = ([0.0, 0.4778, 1.01, 1.458, 1.93, 2.44, 2.93] .± 0.02) ./ 100 
            molar_vol_analysis = MolarVolumeAnalysis(model, pressures_mpa, frac_dilations)
            res_uncertain = molar_vol_analysis.partial_molar_volumes_cm3_mol[3].val

            # with uncertainty and compressibility
            molar_vol_analysis = MolarVolumeAnalysis(model, pressures_mpa, frac_dilations, 1e-5)
            res_uncertain_compressible = molar_vol_analysis.partial_molar_volumes_cm3_mol[3].val

            @test res_no_uncertain ≈ res_uncertain
            @test res_no_uncertain != res_uncertain_compressible

            MolarVolumeAnalysis(DualModeModel(1±1, 4, 5), [1, 2, 3, 4, 5 ± 0.1], [1, 2, 3, 4, 5] ./ 100 .± 0.001, 1e-5)
            
            co2_025_dil_p = [0, 27.3, 66.2, 111.3, 134, 202.6, 256.8, 296.8, 358.2, 407.4, 453.4, 501.4] .* 0.00689476
            co2_025_frac_dil = [0 ± 0, 0.010053518 ± 0.000577937, 0.01908272 ± 0.001100238, 0.02783634 ± 0.001609498, 0.030801961 ± 0.00178267, 0.042389159 ± 0.002462366, 0.050273318 ± 0.002927624, 0.05573762 ± 0.003251388, 0.065036662 ± 0.003804798, 0.070552007 ± 0.004134467, 0.076878536 ± 0.004513929, 0.082775442 ± 0.004868871, ]
            co2_025_mva = MolarVolumeAnalysis(DualModeModel(74.09 ± 0.43, 4.616 ± 0.06, 13.61 ± 0.17), co2_025_dil_p, co2_025_frac_dil, 0)
            co2_025_mva.dfracional_dilation_dp
            co2_025_mva.continuous_dilations


            # no uncertainty, dualmode
            model = DualModeModel(0, 0, 7.8037/(16.698 * 0.101325)) # cc/mpa
            pressures_mpa = [0.0, 4.34, 9.04, 12.8, 17.0, 21.16, 25.15] .* 0.101325
            frac_dilations = [0.0, 0.4778, 1.01, 1.458, 1.93, 2.44, 2.93] ./ 100
            molar_vol_analysis = MolarVolumeAnalysis(model, pressures_mpa, frac_dilations; modeltype=DualModeDilation())
            res_no_uncertain = molar_vol_analysis.partial_molar_volumes_cm3_mol[3]
    
            # with uncertainty, dualmode
            model = DualModeModel(0, 0, 7.8037/(16.698 * 0.101325) ± 0.05) # cc/mpa
            pressures_mpa = ([0.0, 4.34, 9.04, 12.8, 17.0, 21.16, 25.15] .± 0.01) .* 0.101325 
            frac_dilations = ([0.0, 0.4778, 1.01, 1.458, 1.93, 2.44, 2.93] .± 0.02) ./ 100 
            molar_vol_analysis = MolarVolumeAnalysis(model, pressures_mpa, frac_dilations; uncertainty_method = :JackKnife, modeltype=DualModeDilation())
            res_uncertain = molar_vol_analysis.partial_molar_volumes_cm3_mol[3].val
            @test res_no_uncertain == res_uncertain

        end
    end

    # writers
    @testset "Analysis Writer Methods" begin
        results_folder = joinpath(@__DIR__, "writer_output")
        
        if !isdir(results_folder)
            mkdir(results_folder)
        end

        # VHDM 
        path = joinpath(results_folder, "Vant Hoff Dual Mode.xlsx")
        rm(path; force=true)
        write_analysis(VantHoffDualModeAnalysis(isotherms), path)

        # Isosteric Heat
        path = joinpath(results_folder, "Isosteric Heats.xlsx")
        eos_z(p, t) = MembraneEOS.compressibility_factor(MembraneEOS.PR("CH4"), p, t) 
        rm(path; force=true)
        write_analysis(IsostericHeatAnalysis(isotherms), path, name = "Ideal Koros")
        write_analysis(IsostericHeatAnalysis(isotherms, eos_z), path, name = "PREoS Koros")
        
        # Webb Isosteric Heat
        write_analysis(SorptionModels.WebbIsostericHeatAnalysis(isotherms, eos_z), path, name = "Ideal Webb")
        write_analysis(SorptionModels.WebbIsostericHeatAnalysis(isotherms), path, name="PREoS Webb")
        
        
        # Thermo and Mobility Factor 
        path = joinpath(results_folder, "Thermo and Mobility Factors.xlsx")
        rm(path; force=true)
        
        function penetrant_activity(
                p_mpa::Number; 
                T_ref = 273.15, 
                model=Clapeyron.PR(["Carbon Dioxide"])
            )

            Saturation_p_ref = saturation_pressure(model, T_ref)[1]
            mu_ref = Clapeyron.chemical_potential(model, Saturation_p_ref, T_ref)[1]
            p = p_mpa * 1e6
            mu2 = Clapeyron.chemical_potential(model, p, 308.15)[1]
            a = exp((mu2-mu_ref)/(8.314*308.15))
            return a
        end
        
        function CO2_activity(p_MPa::Number)
            return penetrant_activity(p_MPa;T_ref = 273.15, model = Clapeyron.PR(["Carbon Dioxide"]))
        end

        function CH4_activity(p_MPa::Number)
            return penetrant_activity(p_MPa;T_ref = 190, model = Clapeyron.PR(["Methane"]))
        end


        CPIM_CO2_ISOTHERM = IsothermData(
            partial_pressures_mpa=[0.240184901, 0.501416337, 0.866695667, 1.329706241, 1.827098882, 2.444419399, 3.264830659],
            concentrations_cc = [25.11294301, 37.01179749, 49.13345434, 62.09369168, 73.92173216, 85.89637635, 102.769912],
            rho_pol_g_cm3 = 1.285,
            pen_mws_g_mol = 44.009)
        
        CPIMTT_CO2_ISOTHERM = IsothermData(
            partial_pressures_mpa=[0.089419,0.208500733,0.381096791,0.66765144,0.9988614,1.323492497,1.712143839,2.194768515,2.920514824],
            concentrations_cc = [19.89688159, 32.61725119, 44.14091129, 57.54009015, 68.05782353, 78.86482977, 87.59592028, 98.03122486, 118.3535831],
            rho_pol_g_cm3 = 1.244,
            pen_mws_g_mol = 44.009)


        CPIMPPN30_CO2_ISOTHERM = IsothermData(
            partial_pressures_mpa=[0.118208205, 0.331927, 0.566462925, 0.960060366, 1.553171425, 2.310882721, 3.088287685],
            concentrations_cc = [32.11227008, 52.99161609, 67.14573414, 86.1474714, 105.8245687, 126.4041892, 145.6850161],
            rho_pol_g_cm3 = 1.267,
            pen_mws_g_mol = 44.009)


        CPIMPPN30TT_CO2_ISOTHERM = IsothermData(
            partial_pressures_mpa=[0.120226367, 0.244413063, 0.516538011, 0.832279192, 1.171062322, 1.562821945, 2.264679544, 2.933247504],
            concentrations_cc = [36.97127597, 51.61972017, 69.91224668, 83.88229439, 95.89679652, 106.2911621, 124.0630625, 141.289821],
            rho_pol_g_cm3 = 1.244,
            pen_mws_g_mol = 44.009)

        CPIM_CH4_ISOTHERM = IsothermData(
            partial_pressures_mpa=[0.309802563, 0.710304653, 1.222477562, 1.832018444, 2.840712537],
            concentrations_cc = [8.147067101, 14.38986172, 20.52870123, 25.0722395, 32.85525989],
            rho_pol_g_cm3 = 1.285,
            pen_mws_g_mol = 16.043)
        
        CPIMTT_CH4_ISOTHERM = IsothermData(
            partial_pressures_mpa=[0.198645721, 0.490266238, 0.850357598, 1.254338173, 1.855273956, 2.797127243],
            concentrations_cc = [8.012342679, 15.46066477, 21.34320201, 26.24145938, 32.14977817, 40.19937276],
            rho_pol_g_cm3 = 1.244,
            pen_mws_g_mol = 16.043)


        CPIMPPN30_CH4_ISOTHERM = IsothermData(
            partial_pressures_mpa=[0.198068441, 0.483778326, 0.843047932, 1.267734027, 1.772969874, 2.781367538],
            concentrations_cc = [12.11826566, 22.30897037, 30.19814698, 37.13549362, 43.8071852, 53.25404289],
            rho_pol_g_cm3 = 1.267,
            pen_mws_g_mol = 16.043)


        CPIMPPN30TT_CH4_ISOTHERM = IsothermData(
            partial_pressures_mpa=[0.253240639, 0.526341973, 0.849600237, 1.231708745, 1.737795787, 2.592937412],
            concentrations_cc = [17.1136849, 26.44644374, 33.71549609, 40.0423849, 46.45976185, 55.27945851],
            rho_pol_g_cm3 = 1.244,
            pen_mws_g_mol = 16.043)

        
        CPIMs_CO2_Isotherms_Names = ["CPIM_CO2_ISOTHERM", "CPIMTT_CO2_ISOTHERM","CPIMPPN30_CO2_ISOTHERM","CPIMPPN30TT_CO2_ISOTHERM"]
        CPIMs_CH4_Isotherms_Names = ["CPIM_CH4_ISOTHERM","CPIMTT_CH4_ISOTHERM","CPIMPPN30_CH4_ISOTHERM","CPIMPPN30TT_CH4_ISOTHERM"]

        CPIMs_CO2_Isotherms = [CPIM_CO2_ISOTHERM, CPIMTT_CO2_ISOTHERM,CPIMPPN30_CO2_ISOTHERM,CPIMPPN30TT_CO2_ISOTHERM]
        CPIMs_CH4_Isotherms = [CPIM_CH4_ISOTHERM,CPIMTT_CH4_ISOTHERM,CPIMPPN30_CH4_ISOTHERM,CPIMPPN30TT_CH4_ISOTHERM]

        model_types = [DualMode() for _ in 1:length(CPIMs_CO2_Isotherms)]

        CPIM_CH4_DMs = fit_model.(model_types, CPIMs_CH4_Isotherms; uncertainty_method = :JackKnife)
        CPIM_CO2_DMs = fit_model.(model_types, CPIMs_CO2_Isotherms; uncertainty_method = :JackKnife)

        CPIM_CO2_DM = fit_model(DualMode(), CPIM_CO2_ISOTHERM; uncertainty_method = :JackKnife)
        
        
        CPIM_CO2_TFA_CONTINUOUS = ThermodynamicFactorAnalysis(CPIM_CO2_ISOTHERM, CPIM_CO2_DM, x -> CO2_activity(x))
        CPIMs_CO2_TFA_CONTINUOUS = ThermodynamicFactorAnalysis.(CPIMs_CO2_Isotherms, CPIM_CO2_DMs,x->CO2_activity(x))
        CPIMs_CH4_TFA_CONTINUOUS = ThermodynamicFactorAnalysis.(CPIMs_CH4_Isotherms, CPIM_CH4_DMs, x->CH4_activity(x))
        
        
        [write_analysis(CPIMs_CO2_TFA_CONTINUOUS[i],path, name = CPIMs_CO2_Isotherms_Names[i])  for i in eachindex(CPIMs_CO2_Isotherms_Names,CPIMs_CO2_TFA_CONTINUOUS)]

        [write_analysis(CPIMs_CH4_TFA_CONTINUOUS[i],path, name = CPIMs_CH4_Isotherms_Names[i])  for i in eachindex(CPIMs_CH4_Isotherms_Names,CPIMs_CH4_TFA_CONTINUOUS)]


        # begin diffusivities
        CPIM_CO2_Diffusivities = [2.51E-07 ± 1.56E-08, 2.64E-07 ± 1.54E-08, 2.87E-07 ± 1.54E-08, 3.21E-07 ± 1.56E-08, 3.35E-07 ± 1.57E-08, 4.33E-07 ± 1.66E-08, 6.06E-07 ± 1.97E-08, 7.57E-07 ± 2.31E-08, 9.18E-07 ± 2.70E-08, 1.11E-06 ± 3.19E-08, 1.23E-06 ± 3.51E-08, 1.35E-06 ± 3.82E-08]
        CPIM_CO2_D_Pressures = [0.02, 0.04, 0.07, 0.11, 0.13, 0.29, 0.63, 1.01, 1.46, 2.18, 2.55, 3.01]

        CPIMTT_CO2_Diffusivities =[1.23E-07, 1.62E-07, 2.42E-07, 2.92E-07, 3.99E-07, 4.86E-07, 5.39E-07, 6.01E-07, 6.61E-07, 6.92E-07, 7.40E-07, 7.56E-07, 7.95E-07]
        CPIMTT_CO2_D_Pressures = [0.026108491, 0.096498915, 0.329691673, 0.51115978, 1.013709993, 1.533341029, 2.053176956, 2.474938589, 3.074988743, 3.544684124, 4.049162551, 4.528244384, 5.018626422]

        CPIMPPN30_CO2_Diffusivities = [5.70E-08, 6.28E-08, 7.51E-08]
        CPIMPPN30_CO2_D_Pressures = [0.10575012, 0.149904497, 0.236765912]

        CPIMPPN30TT_CO2_Diffusivities = [2.56E-07, 2.67E-07, 2.84E-07, 3.15E-07, 3.35E-07, 3.70E-07, 5.83E-07, 7.91E-07, 9.96E-07, 1.17E-06, 1.31E-06, 1.36E-06, 1.45E-06, 1.57E-06, 1.67E-06, 1.70E-06]
        CPIMPPN30TT_CO2_D_Pressures = [0.023958397, 0.051496768, 0.071641224, 0.098420661, 0.122264797, 0.171305182, 0.493184625, 0.897479829, 1.488763245, 2.030096644, 2.46419457, 2.96169734, 3.431452286, 4.017734468, 4.553939642, 5.012233393]

        CPIM_CH4_Diffusivities = [6.00E-08, 7.21E-08, 8.27E-08, 9.68E-08]
        CPIM_CH4_D_Pressures = [0.108381747, 0.352283335, 0.632820204, 0.995498696]

        CPIMTT_CH4_Diffusivities =[0.106659394, 0.348679577, 0.681648981, 1.020196443, 1.520182062]
        CPIMTT_CH4_D_Pressures = [3.65E-08, 5.14E-08, 6.19E-08, 6.97E-08, 8.18E-08]

        CPIMPPN30_CH4_Diffusivities = [1.00E-08, 1.40E-08, 1.61E-08, 1.84E-08]
        CPIMPPN30_CH4_D_Pressures = [0.113271481, 0.522281552, 0.739791262, 1.040099094]

        CPIMPPN30TT_CH4_Diffusivities = [5.49E-08, 6.94E-08, 9.04E-08, 1.09E-07]
        CPIMPPN30TT_CH4_D_Pressures = [0.104457033, 0.317433866, 0.609517268, 0.965470726]

        CPIMS_CO2_D_PRESSURES = [CPIM_CO2_D_Pressures,CPIMTT_CO2_D_Pressures,CPIMPPN30_CO2_D_Pressures,CPIMPPN30TT_CO2_D_Pressures]
        CPIMS_CO2_DIFFUSIVITIES = [CPIM_CO2_Diffusivities,CPIMTT_CO2_Diffusivities,CPIMPPN30_CO2_Diffusivities,CPIMPPN30TT_CO2_Diffusivities]
        
        CPIM_CO2_L_ANALYSES = [MobilityFactorAnalysis(
            CPIMS_CO2_DIFFUSIVITIES[i],
            CPIMS_CO2_D_PRESSURES[i],
            polymer_density(CPIMs_CO2_Isotherms[i]),
            44.009,
            CPIM_CO2_DMs[i], CO2_activity)
            for i in eachindex(CPIMs_CO2_Isotherms_Names)]
        
        [write_analysis(CPIM_CO2_L_ANALYSES[i],path, name = CPIMs_CO2_Isotherms_Names[i])  for i in eachindex(CPIMs_CO2_Isotherms_Names)]


        CPIM_CO2_L_ANALYSIS = MobilityFactorAnalysis(CPIM_CO2_Diffusivities, CPIM_CO2_D_Pressures, 1.285, 44.009, CPIM_CO2_DM, penetrant_activity)
      
        write_analysis(CPIM_CO2_L_ANALYSIS, path, name="CPIM_CO2_L_ANALYSIS")
        
        [write_analysis(CPIMs_CO2_TFA_CONTINUOUS[i],path, name = CPIMs_CO2_Isotherms_Names[i])  for i in eachindex(CPIMs_CO2_Isotherms_Names,CPIMs_CO2_TFA_CONTINUOUS)]

        

        CPIM_CO2_T_ANALYSIS = ThermodynamicFactorAnalysis(
            collect(
                range(
                    0.001, 
                    5, 
                    step=0.02
                )
            ), 
            1.285, 44.009, CPIM_CO2_DM, CO2_activity)


        write_analysis(CPIM_CO2_T_ANALYSIS, path, name="CPIM_CO2_T_ANALYSIS")
       
        CPIM_CO2_T_ANALYSES = [ThermodynamicFactorAnalysis(
            collect(
                range(
                    0.001, 
                    5, 
                    step=0.02
                )
            ), 
            polymer_density(CPIMs_CO2_Isotherms[i]), 44.009, CPIM_CO2_DMs[i], CO2_activity) for i in eachindex(CPIM_CO2_DMs)]

        [write_analysis(CPIM_CO2_T_ANALYSES[i],path, name = CPIMs_CO2_Isotherms_Names[i])  for i in eachindex(CPIMs_CO2_Isotherms_Names,CPIM_CO2_T_ANALYSES)]

        # Partial Immobilization 
        path = joinpath(results_folder, "Partial Immobilization Model.xlsx")
        rm(path; force=true)
        model = DualModeModel(56.8, 7.4, 26.1; use_fugacity=true)
        permeabilities = [1221, 1091, 1038]
        pressures_mpa = [0.298, 0.637, 0.974]
        result = PartialImmobilizationModel(model, pressures_mpa, permeabilities)
        write_analysis(result, path)

        # Zimm-Lundberg
        path = joinpath(results_folder, "Zimm-Lundberg Analysis.xlsx")
        rm(path; force=true)
            # water isotherm taken from from TPBO-0.25 at 25C
        activities=[0.068087235, 0.129051621, 0.185054022, 0.226470588, 0.26212485, 0.317987195, 0.417687075, 0.54869948, 0.641314436]
        concentrations_cc = [4.553554612, 7.361278518, 9.973016055, 12.76538706, 14.32931957, 17.19635897, 22.19538967, 28.84311218, 33.56853954]
        activities = activities .± (activities .* 0.001); concentrations_cc = concentrations_cc .± (concentrations_cc .* 0.001)
        mol_vol = 18 ± 0.2 # cm3/mol
        isotherm = IsothermData(; activities, concentrations_cc)
        gabmodel = fit_model(GAB(), isotherm)
        result = ZimmLundbergAnalysis(gabmodel, activities, mol_vol)
        write_analysis(result, path)

        # Dual Mode Desorption

        isotherm_desorption_with_errs = IsothermData(; 
            partial_pressures_mpa = [0.241333352, 0.600763584, 1.04806673, 1.466095481, 1.951571285, 2.499847618, 3.142683031, 2.60974313, 1.199714642, 0.575992209, 0.30991402, 0.145032338] .± 0.01,
            concentrations_cc = [41.48924079, 62.79313671, 77.9590348, 88.08019013, 96.61909374, 105.7659302, 114.9523482, 111.9833899, 98.93545196, 81.21343867, 66.62816092, 50.92125421] .± 0.01
        )
        path = joinpath(results_folder, "Dual Mode Desorption.xlsx")
        rm(path; force=true)
        dmda = DualModeDesorption(isotherm_desorption_with_errs; uncertainty_method=:JackKnife)
        write_analysis(dmda, path)

        # Partial Molar Volumes
        path = joinpath(results_folder, "Molar Volume Analysis.xlsx")
        rm(path; force=true)
        model = DualModeModel(0, 0, 7.8037/(16.698 * 0.101325)) # cc/mpa
        pressures_mpa = [0.0, 4.34, 9.04, 12.8, 17.0, 21.16, 25.15] .* 0.101325
        frac_dilations = [0.0, 0.4778, 1.01, 1.458, 1.93, 2.44, 2.93] ./ 100
        molar_vol_analysis = MolarVolumeAnalysis(model, pressures_mpa, frac_dilations; uncertainty_method=nothing)
        write_analysis(molar_vol_analysis, path)
       
        model = DualModeModel(0, 0, 7.8037/(16.698 * 0.101325) ± 0.05) # cc/mpa
        pressures_mpa = ([0.0, 4.34, 9.04, 12.8, 17.0, 21.16, 25.15] .± 0.01) .* 0.101325 
        frac_dilations = ([0.0, 0.4778, 1.01, 1.458, 1.93, 2.44, 2.93] .± 0.02) ./ 100 
        molar_vol_analysis = MolarVolumeAnalysis(model, pressures_mpa, frac_dilations)
        write_analysis(molar_vol_analysis, path; name="With hessian an uncertain input")

        model = DualModeModel(0, 0, 7.8037/(16.698 * 0.101325) ± 0.05) # cc/mpa
        pressures_mpa = ([0.0, 4.34, 9.04, 12.8, 17.0, 21.16, 25.15] .± 0.01) .* 0.101325 
        frac_dilations = ([0.0, 0.4778, 1.01, 1.458, 1.93, 2.44, 2.93] .± 0.02) ./ 100 
        molar_vol_analysis = MolarVolumeAnalysis(model, pressures_mpa, frac_dilations, uncertainty_method=:JackKnife)
        write_analysis(molar_vol_analysis, path; name="Also with jackknifing!")

            # one particular error prone case        
        co2_000_p = [0.241333352, 0.600763584, 1.04806673, 1.466095481, 1.951571285, 2.499847618, 3.142683031]
        co2_000_c = [41.48924079, 62.79313671, 77.9590348, 88.08019013, 96.61909374, 105.7659302, 114.9523482]
        co2_000_dil_p = [0, 44.5, 104.6, 153.4, 231, 336.8, 381.5, 433.5, 470.3] .* 0.00689476
        co2_000_frac_dil = [0 ± 0, 0.012127622 ± 0.00044235, 0.026889891 ± 0.000985522, 0.036407288 ± 0.001338418, 0.048826916 ± 0.001802068, 0.06476425 ± 0.00240217, 0.071486539 ± 0.002656993, 0.078237062 ± 0.00291389, 0.08335314 ± 0.003109253]
        iso = IsothermData(partial_pressures_mpa = co2_000_p, concentrations_cc = co2_000_c)
        dmmodel = fit_model(DualMode(), iso; uncertainty_method = :JackKnife)
        co2_000_mva_many_param = MolarVolumeAnalysis(dmmodel, co2_000_dil_p, co2_000_frac_dil, 0; uncertainty_method=:JackKnife, n_params=4)
        co2_000_mva_few_params = MolarVolumeAnalysis(dmmodel, co2_000_dil_p, co2_000_frac_dil, 0; uncertainty_method=:JackKnife, n_params=3)
        @test Measurements.uncertainty.(co2_000_mva_few_params.continuous_dilations) != Measurements.uncertainty.(co2_000_mva_many_param.continuous_dilations) 
        @test Measurements.uncertainty.(co2_000_mva_many_param.continuous_dilations) != Measurements.uncertainty.(co2_000_mva_few_params.continuous_dilations) 

        molar_vol_analysis = MolarVolumeAnalysis(model, pressures_mpa, frac_dilations, uncertainty_method=:JackKnife, modeltype=DualModeDilation())
        write_analysis(molar_vol_analysis, path; name="dualmode dilation & jackknife")
        

    end

    @testset "Diagnostic Methods" begin
        results_folder = joinpath(@__DIR__, "diagnostics_output")
        
        # NELF
        isotherms = [tpbo_ch4_5c, tpbo_ch4_20c, tpbo_ch4_35c, tpbo_co2_5c, tpbo_co2_20c, tpbo_co2_35c, tpbo_co2_50c, tpbo_n2_5c, tpbo_n2_50c]
        char_co2 = [630, 300, 1.515, 44]
        char_ch4 = [250, 215, 0.500, 16.04]
        char_n2 = [160, 145, 0.943, 28.01]
        bulk_phase_char_params = [char_ch4, char_ch4, char_ch4, char_co2, char_co2, char_co2, char_co2, char_n2, char_n2]
        error_plot = nelf_characteristic_parameter_error_map(isotherms, bulk_phase_char_params, verbose=false)
        for val in error_plot
            if isnan(val)
                @show "Something very bad happened in the NELF diagnostics"
            end
        end
        # savefig(error_plot, joinpath(results_folder, "TPBO_25 error map.png"))
    
        # DGRPT
        polymer = "PC"
        penetrant = "CO2"
        kij = [0 -0.007; -0.007 0]
        bulk_phase_eos = MembraneEOS.SL([penetrant])
        polymer_phase_eos = MembraneEOS.SL([polymer, penetrant], kij)
        density = 1.197850471   #g/cm3    
        temperature = 308.15
        dgrptmodel = DGRPTModel(bulk_phase_eos, polymer_phase_eos, density)
            
        ps = [0.1, 1, 2, 3, 3.5]

        dgrpt_concs_pure_co2_1_term = [predict_concentration(dgrptmodel, temperature, p, [1.0]; taylor_series_order=1, units=:frac)[2] for p in ps]
        dgrpt_1_term_mass_fracs = [[1-mfrac, mfrac] for mfrac in dgrpt_concs_pure_co2_1_term]
        dgrpt_dens_pure_co2_1_term = [polymer_density(dgrptmodel, temperature, p, [1]; taylor_series_order=1) for p in ps]
        
        function plot_density_target_function!(plt, model, temp, polphase_mass_fracs)
            target_roots = SorptionModels.make_roots_polymer_density_target(model, temp, polphase_mass_fracs; taylor_series_order=1)
            max_polymer_density = SorptionModels.density_upper_bound(model.polymer_model, polphase_mass_fracs) * polphase_mass_fracs[1]
            tested_densities = range(eps(), max_polymer_density, 400)
            plot!(plt, tested_densities, target_roots.(tested_densities), label="pol mfrac: " * string(polphase_mass_fracs[1]))
        end

        target_easy_roots = SorptionModels.make_roots_polymer_density_target(dgrptmodel, temperature, dgrpt_1_term_mass_fracs[1]; taylor_series_order=1)
        target_easy_optim = SorptionModels.make_optim_polymer_density_target(dgrptmodel, temperature, dgrpt_1_term_mass_fracs[1]; taylor_series_order=1)
        max_polymer_density_easy = SorptionModels.density_upper_bound(dgrptmodel.polymer_model, dgrpt_1_term_mass_fracs[1]) * dgrpt_1_term_mass_fracs[1][1]
        tested_densities = range(eps(), max_polymer_density_easy, 400)
        resplot = plot(tested_densities, target_easy_roots.(tested_densities), label="roots")
        # plot!(resplot, tested_densities, target_easy_optim.(tested_densities), label="optim")
        resplot = plot()
        # plot_density_target_function!(resplot, dgrptmodel, temperature, dgrpt_1_term_mass_fracs[2])
        for mfracs in dgrpt_1_term_mass_fracs
            plot_density_target_function!(resplot, dgrptmodel, temperature, mfracs)
        end
        savefig(resplot, joinpath(results_folder, "DGRPT PC-CO2 densitiy targets.png"))


        # dgrpt_concs_pure_co2_2_terms = [predict_concentration(dgrptmodel, temperature, p, [1.0]; taylor_series_order=2)[1] for p in ps]
        # dgrpt_concs_pure_co2_5_terms = [predict_concentration(dgrptmodel, temperature, p, [1.0]; taylor_series_order=5)[1] for p in ps]
        # @show dgrpt_concs_pure_co2_1_term
        # @show dgrpt_dens_pure_co2_1_term
        
    end

    @testset "Misc Compatibility" begin
        # Polynomials.jl and Measurements.jl should play nice together. JuliaPhysics/Measurements.jl#134 caused 1.9 to be incompatible with Polynomials.jl
        #   one(x) in measurements returns the underlying measurement type rather than `measurement`, so until polynomials fixes https://github.com/JuliaMath/Polynomials.jl/issues/501, measurements must be restricted to pre-1.9 versions.
        x = [0, 27.3, 66.2, 111.3, 134, 202.6, 256.8, 296.8, 358.2, 407.4, 453.4, 501.4] .* 0.00689476
        y = [0 ± 0, 0.010053518 ± 0.000577937, 0.01908272 ± 0.001100238, 0.02783634 ± 0.001609498, 0.030801961 ± 0.00178267, 0.042389159 ± 0.002462366, 0.050273318 ± 0.002927624, 0.05573762 ± 0.003251388, 0.065036662 ± 0.003804798, 0.070552007 ± 0.004134467, 0.076878536 ± 0.004513929, 0.082775442 ± 0.004868871, ]
        SorptionModels.Polynomials.fit(x, y, 4)

    end
end
nothing
