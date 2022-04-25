using Revise
using SorptionModels
using Test
using Measurements
using MembraneBase

precision = 5

@testset "SorptionModels.jl" begin

    
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
        dmfittings = fit_dualmode_model(isofit)
        
        @test round(dmfittings.ch; digits=precision) == round(6.807641216606124; digits=precision) &&
            round(dmfittings.b; digits=precision) == round(3.376206925261848; digits=precision) && 
            round(dmfittings.kd; digits=precision) == round(5.556359995999332; digits=precision)


        # test jackknife method (same values as the first dmfittings)
        dmfittings_3 = fit_dualmode_model(isofit; uncertainty_method=:JackKnife)
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
        @test model.ch != model_different.ch

    # GAB
        acts =  [0, 0.023, 0.083, 0.161, 0.202,  0.263, 0.327 ]  # this is real 1-propanol sorption in TPBO-0.25 at 25C
        concs = [0, 2.253, 3.925, 6.126, 8.091, 12.363, 20.261]


        # test basic fitting consistency 
        gmfitting = fit_gab_model(acts, concs)
        @test gmfitting.cp == 4.7515508907549435 && 
            gmfitting.k ==  2.365371628472392 && 
            gmfitting.a == 9.807108219804844
            
        # test checking for DimensionMismatch
        @test_throws DimensionMismatch fit_gab_model([1, 2, 3], [1, 2, 3, 4])

        # test fitting with Measurement types
        acts_2 = [0 ± 0.1, 0.023 ± 0.1, 0.083 ± 0.1, 0.161 ± 0.1, 0.202 ± 0.1, 0.263 ± 0.1, 0.327 ± 0.1]  # this is real 1-propanol sorption in TPBO-0.25 at 25C
        concs_2 = [0 ± 0.5, 2.253 ± 0.5, 3.925 ± 0.5, 6.126 ± 0.5, 8.091 ± 0.5, 12.363 ± 0.5, 20.261 ± 0.5]
        
        gmfitting_2 = fit_gab_model(acts_2, concs_2)
        @test round(gmfitting_2.cp; digits = 5) == round(4.7515508907549435; digits = 5) && 
            round(gmfitting_2.k; digits = 5) ==  round(2.36537162847239; digits = 5) && 
            round(gmfitting_2.a; digits = 5) == round(9.807108219804844; digits = 5)
        
        # test fitting with JackKnife uncertainty_method
        gmfitting_3 = fit_gab_model(acts_2, concs_2; uncertainty_method=:JackKnife)
        @test gmfitting_3.cp.err ≈ 0.49739433319622156 && 
            gmfitting_3.k.err ≈ 0.1236803472661029 && 
            gmfitting_3.a.err ≈ 4.09125368150671
    end

    @testset "Fundamental Sorption Models" begin
        using MembraneEOS
        polymer = "PC"
        penetrant = "CO2"
        kij = [0 -0.007; -0.007 0]
        bulk_phase_eos = SL([penetrant])
        polymer_phase_eos = SL([polymer, penetrant], kij)
        density = 1.197850471   #g/cm3    
        temperature = 308.15
        pressures = [0, 0.18, 0.38, 0.64, 0.94, 1.23, 1.44]
        expected_mass_fracs = [5.51E-06, 0.008294923, 0.014447025, 0.020467468, 0.026066798, 0.030734723, 0.033827052]
        expected_concs_cc_cc = [0, 5.094493596, 8.910071012, 12.66689539, 16.17497812, 19.10613428, 21.05001223]
        
        # NELF    
            ksw = 0.0102            # 1/MPa
            nelfmodel = NELFModel(bulk_phase_eos, polymer_phase_eos, density, ksw)
            nelf_concs_pure_co2 = [predict_concentration(nelfmodel, temperature, p, [1.0])[1] for p in pressures]
        
            # @show [PolymerMembranes.bulk_phase_chemical_potential(nelfmodel, temperature, p, [1])[1] for p in pressures]
        
            penetrants = ["CO2", "O2"]
            kij_ternary = [0      -0.007 0.0 ; 
                           -0.007 0      0.0 ; 
                           0      0.0    0.0 ]
            ksw_ternary = [0.0102, 0.0]            # 1/MPa
            bulk_phase_eos_ternary = SL(penetrants)
            polymer_phase_eos_ternary = SL([polymer, penetrants...], kij_ternary)
            nelfmodel_ternary = NELFModel(bulk_phase_eos_ternary, polymer_phase_eos_ternary, density, ksw_ternary)
            nelf_concs_co2_mix = [predict_concentration(nelfmodel_ternary, temperature, p, [0.5, 0.5])[1] for p in pressures]
            @test nelf_concs_co2_mix[3] != nelf_concs_pure_co2[3]
            
            nelf_concs_co2_psuedo = [predict_concentration(nelfmodel_ternary, temperature, p, [1.0, 0])[1] for p in pressures]
            @test nelf_concs_co2_psuedo[3] ≈ nelf_concs_pure_co2[3]
        # DGRPT

            # dgrptmodel = DGRPTModel(bulk_phase_eos, polymer_phase_eos, density)
            # @show dgrpt_concs_pure_co2 = [predict_concentration(dgrptmodel, temperature, p, [1.0])[1] for p in pressures]

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
        exp_data = TransientStepData(
            [1,       2,       3,        4,       5,      6,      7,      8,      9,      10,    15,    20,    25,    30,    40,   50,    70,   90,     100], 
            [0.00194, 0.00515, 0.009107, 0.01359, 0.0184, 0.0237, 0.0291, 0.0348, 0.0407, 0.046, 0.077, 0.109, 0.141, 0.171, 0.22, 0.278, 0.36, 0.4364, 0.4671]            
        )
        fick_model_fit = fit_transient_sorption_model(exp_data, FickianSorption())
        
        bh_model_fit = fit_transient_sorption_model(exp_data, BerensHopfenbergSorption())
    
        mbh_model_fit = fit_transient_sorption_model(exp_data, ModifiedBerensHopfenbergSorption())
        
        semi_thickness_cm = 0.02 ± 0.001

        fick_d = get_diffusivity(fick_model_fit, semi_thickness_cm)
        @test fick_d.val ≈ 4.779020339107122e-7
        @test fick_d.err ≈ 4.7790203391071216e-8

        mbh_d = get_diffusivity(mbh_model_fit, semi_thickness_cm)
        @test mbh_d.val ≈ 2.251436090950463e-15
        @test mbh_d.err ≈ 2.251436090950463e-16

        # boostrap and uncertainty
        fick_model_fit_with_err = fit_transient_sorption_model(exp_data, FickianSorption(); uncertainty_method=:Bootstrap)
        fick_d_with_err = get_diffusivity(fick_model_fit_with_err, semi_thickness_cm)
        @test fick_d_with_err.val ≈ fick_d.val
        # @test fick_d_with_err.err ≈ 7.950417186869458e-8 # errors are random due to low sample size
        
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
        fugacities_mpa =[0.27337983, 0.684971536, 1.154961289, 1.651558222, 2.128304942, 2.599274758, 2.985937811]
    )
    isotherms = [isotherm_1, isotherm_2, isotherm_3, isotherm_4]
    @testset "Analysis Functions" begin
        
        # VHDM analysis
        @testset "Van't Hoff Dual Mode Fitting" begin    
            vhdm_analysis = VantHoffDualModeAnalysis(isotherms)
            vhdm_analysis_but_with_fugacity = VantHoffDualModeAnalysis(isotherms; use_fugacity=true)
            @test vhdm_analysis_but_with_fugacity.final_models[1].use_fugacity == true
            @test vhdm_analysis_but_with_fugacity.final_models[2].ch.val ≈ 57.13799466676326
            @test vhdm_analysis_but_with_fugacity.final_models[2].ch.err ≈ 29.838724853992588
            @test vhdm_analysis.final_models[2].ch.val ≈ 56.37550904605368
        end

        # Isosteric Heat
        @testset "Isosteric Heat Analysis" begin    
            ish_analysis = IsostericHeatAnalysis(isotherms)
            ish_analysis_but_with_vhdm = IsostericHeatAnalysis(isotherms; use_vant_hoff_constraints=true)

            conc_no_constraint = ish_analysis.sampled_concentrations[4]
            conc_with_constraint = ish_analysis_but_with_vhdm.sampled_concentrations[4]
            @test conc_no_constraint == conc_with_constraint
            @test ish_analysis.isosteric_heat_at_conc[4].val ≈ -18998.600511017383
            @test ish_analysis_but_with_vhdm.isosteric_heat_at_conc[4].val ≈ -18118.53623433742
        end

        # Mobility Factor
        @testset "Mobility Factor" begin
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
            mob_fact_analysis = MobilityFactorAnalysis(iso, dif)
            @test mob_fact_analysis.kinetic_factors[3].val ≈ 8.07131661298753e-7
            @test mob_fact_analysis.thermodynamic_factors[3].val ≈ 1.0273999147371586
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
    end
    # writers
    @testset "Analysis Writer Methods" begin
        results_folder = joinpath(@__DIR__, "writer_output")

        # VHDM 
        path = joinpath(results_folder, "Vant Hoff Dual Mode Analysis.xlsx")
        rm(path; force=true)
        write_analysis(VantHoffDualModeAnalysis(isotherms), path)

        # Isosteric Heat
        path = joinpath(results_folder, "Isosteric Heat Analysis.xlsx")
        rm(path; force=true)
        write_analysis(IsostericHeatAnalysis(isotherms), path)

        # no mobility factor implemented yet
        
        # Partial Immobilization 
        path = joinpath(results_folder, "Partial Immobilization Model.csv")
        rm(path; force=true)
        model = DualModeModel(56.8, 7.4, 26.1; use_fugacity=true)
        permeabilities = [1221, 1091, 1038]
        pressures_mpa = [0.298, 0.637, 0.974]
        result = PartialImmobilizationModel(model, pressures_mpa, permeabilities)
        write_analysis(result, path)

    end
end
nothing



