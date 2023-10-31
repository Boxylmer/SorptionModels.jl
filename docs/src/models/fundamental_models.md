# Fundamental Sorption Models

## NELF

The NELF model is a fundamental model that describes sorption in polymers through lattice fluid equations of state. The actual underlying equation of state 

```@autodocs
Modules = [SorptionModels]
Pages = ["SorptionModels/NELF.jl"]
```

## DGRPT

The dry glass reference perturbation theory (DGRPT) model attempts to improve the NELF model by predicting polymer swelling rather than being supplied (or fitting) a swelling coefficient for each component in the polymer phase. 

```@autodocs
Modules = [SorptionModels]
Pages = ["SorptionModels/DGRPT.jl"]
```
