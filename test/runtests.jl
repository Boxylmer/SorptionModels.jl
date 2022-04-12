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
        isofit = IsothermData(; partial_pressures_mpa = [0, 0.03, 0.15, 0.6, 0.9, 1.2, 1.5], concentrations_cc = [0, 1, 3, 8, 10, 12.2, 14])
        dmfittings = fit_dualmode_model(isofit)
        
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

        # dmfittings_3 = DualMode.fit_model_to_isotherm(isofit; uncertainty_method=:Bootstrap)
        # @test dmfittings_3.ch.err == 0.756291283433861 &&
        # dmfittings_3.b.err == 1.0368748268249859 && 
        # dmfittings_3.kd.err == 0.3217790802072902


        model = SorptionModels.DualModeModel(2.3, 4.5, 0.1)
        pressure = 5.2
        concentration = predict_concentration(model, pressure)
        recalculated_pressure = SorptionModels.predict_pressure(model, concentration)
        @test pressure == round(recalculated_pressure; digits=precision)

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
        @test gmfitting_2.cp == 4.75152770782204 && 
            gmfitting_2.k ==  2.365375206030592 && 
            gmfitting_2.a == 9.807283609034045

        # test fitting with JackKnife uncertainty_method
        gmfitting_3 = fit_gab_model(acts_2, concs_2; uncertainty_method=:JackKnife)
        @test gmfitting_3.cp.err == 0.4973699474088062 && 
            gmfitting_3.k.err ==  0.12367885393402438 && 
            gmfitting_3.a.err == 4.090851660381514
        
    end

    @testset "Fundamental Sorption Models" begin
        using MembraneEOS
        # NELF
            polymer = "PC"
            penetrant = "CO2"
            kij = [0 0; 0 0]
            ksw = 0.0102            # 1/MPa
            density = 1.197850471   #g/cm3    
            bulk_phase_eos = SL([penetrant])
            polymer_phase_eos = SL([polymer, penetrant], kij)
            # nelfmodel = NELFModel(bulk_phase_eos, polymer_phase_eos, density, ksw)
        
            temperature = 308.15
            pressures = [1e-13, 0.18, 0.38, 0.64, 0.94, 1.23, 1.44]
            expected_mass_fracs = [5.51E-06, 0.008294923, 0.014447025, 0.020467468, 0.026066798, 0.030734723, 0.033827052]
            expected_concs_cc_cc = [0.003361991, 5.094493596, 8.910071012, 12.66689539, 16.17497812, 19.10613428, 21.05001223]
            expected_co2_polymer_phase_μ = [-51.54243958, -32.12143859, -30.22735524, -28.91840386, -27.96456057, -27.30597109, -26.92427275] * 1000
            expected_co2_bulk_phase_μ = [-57.20578109, -32.1214365, -30.22735377, -28.91840261, -27.96455943, -27.30596998, -26.92427169] * 1000
            # acquired_μ = [predict_concentration(nelfmodel, temperature, p, [1.0])[1] for p in pressures]
            # acquired_conc = Vector{Any}(undef, length(expected_mass_fracs))
            # for idx in eachindex(expected_mass_fracs, pressures)
            #     pressure = pressures[idx]
            #     polymer_mass_fraction = 1 - sum([expected_mass_fracs[idx]])
            #     polymer_phase_mass_fractions = vcat(polymer_mass_fraction, [expected_mass_fracs[idx]])
                
            #     # acquired_μ[idx] = PolymerMembranes.polymer_phase_activities(nelfmodel, temperature, pressure, [1], polymer_phase_mass_fractions)[2]
            #     acquired_conc[idx] = predict_concentration(nelfmodel, temperature, pressure, [1], 
            #         [PolymerMembranes.molecular_weight(penetrant)])
            # end
            # @show expected_concs_cc_cc
            # @show acquired_conc
        
            # @show expected_co2_bulk_phase_μ
            # @show [PolymerMembranes.bulk_phase_chemical_potential(nelfmodel, temperature, p, [1])[1] for p in pressures]
        
            # polymer= ChemicalParameters("PC")
            # penetrants = [ChemicalParameters("CO2"), ChemicalParameters("O2")]
            # kij = [0     0     0.0  ; 
            #        0     0     -1.0; 
            #        0     -1. 0    ]
            # ksw = [0.0, 0]            # 1/MPa
            # density = 1.2   #g/cm3
            # nelfmodel = NELFModel(penetrants, polymer, kij, density, ksw)
            # pressures = [1.00E-05, 0.086206897, 0.172413793, 0.25862069, 0.344827586, 0.431034483, 0.517241379, 0.603448276, 0.689655172, 0.775862069, 0.862068966]
            # expected_o2_μ = [-51.40654702, -28.19935774, -26.43327277, -25.40425664, -24.67700853, -24.11510845, -23.65778994, -23.27263698, -22.94029887, -22.64829267, -22.38809381] * 1000
            # expected_co2_μ = [-57.20648535, -34.01571204, -32.26627046, -31.25413255, -30.54400748, -29.99948591, -29.55981291, -29.19258481, -28.87846425, -28.60498276, -28.36363139] * 1000
            # results = [PolymerMembranes.bulk_phase_chemical_potential(nelfmodel, temperature, p * 2, [0.5, 0.5]) for p in pressures]
            # acquired_o2_μ = [r[2] for r in results]
            # acquired_co2_μ = [r[1] for r in results]
            # @show expected_o2_μ
            # @show acquired_o2_μ
        
            # @show expected_co2_μ
            # @show acquired_co2_μ
            # @show abs.(expected_concs_cc_cc .- acquired_concs_cc_cc) ./ expected_concs_cc_cc
            # @show [predict_concentration(nelfmodel, 308.15, p, [1.0]) for p in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]*0.101325]
    end

    @testset "Transient Sorption Models" begin
        
        fickian_model = FickianSorptionModel(1, 0.01)
        bh_model = BerensHopfenbergSorptionModel(0.7, 0.01, 0.3, 0.001)
        mbh_model = ModifiedBerensHopfenbergSorptionModel(0.7, 0.01, 0.3, 0.001, 0.05)
        fickian_prediction = predict_sorption(fickian_model, 10)
        bh_prediction = predict_sorption(bh_model, 10)
        mbh_prediction = predict_sorption(mbh_model, 10)
    
        # println(fickian_prediction)
        # println(bh_prediction)
        # println(mbh_prediction)
    
        time_vector = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30, 40, 50, 70, 90, 100]
        fick_prediction_vector = predict_sorption(fickian_model, time_vector)
        bh_prediction_vector = predict_sorption(bh_model, time_vector)
        mbh_prediction_vector = predict_sorption(mbh_model, time_vector)
        
        
        # testing fitting without errors (no true tests)
        exp_data = TransientStepData(
            [1,       2,       3,        4,       5,      6,      7,      8,      9,      10,    15,    20,    25,    30,    40,   50,    70,   90,     100], 
            [0.00194, 0.00515, 0.009107, 0.01359, 0.0184, 0.0237, 0.0291, 0.0348, 0.0407, 0.046, 0.077, 0.109, 0.141, 0.171, 0.22, 0.278, 0.36, 0.4364, 0.4671]            
        )
        fick_model_fit = fit_transient_sorption_model(exp_data, :FickianSorptionModel)
    
        bh_model_fit = fit_transient_sorption_model(exp_data, :BerensHopfenbergSorptionModel)
    
        mbh_model_fit = fit_transient_sorption_model(exp_data, :ModifiedBerensHopfenbergSorptionModel)
        

    end
end




