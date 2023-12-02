
v★(P★, T★,) = 8.31446261815324 * T★ / P★ / 1000000 # J / (mol*K) * K / mpa -> pa * m3 / (mol * mpa) ->  need to divide by 1000000 to get m3/mol
ϵ★(T★) = 8.31446261815324 * T★ # J / (mol * K) * K -> J/mol 
r(P★, T★, ρ★, mw) = mw * (P★ * 1000000) / (R * T★ * (ρ★ / 1e-6)) # g/mol * mpa * 1000000 pa/mpa / ((j/mol*K) * K * g/(cm3 / 1e-6 m3/cm3)) -> unitless