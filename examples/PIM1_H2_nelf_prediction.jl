using SorptionModels
using Clapeyron



pim1_p★ = 524	
pim1_t★ = 872	
pim1_ρ★ = 1.438	
pim1_dry_ρ = 1.143

h2_p★ = 37
h2_t★ = 46
h2_ρ★ = 0.078


pol_phase_eos = SorptionModels.construct_binary_sl_eosmodel([pim1_p★, pim1_t★, pim1_ρ★, 1e20], [h2_p★, h2_t★, h2_ρ★, 2.016], 0)

nelfmodel = NELFModel(pol_phase_eos, pim1_dry_ρ)

pred_ps_mpa = collect(range(0, 2, length=20))
pred_cs_cc = [predict_concentration(nelfmodel, 308.15, p)[1] for p in pred_ps_mpa]




