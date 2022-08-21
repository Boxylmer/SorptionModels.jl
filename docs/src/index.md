# SorptionModels.jl

SorptionModels.jl provides a number of models which describe solubility in polymers. They are broken down into three categories:

### Empirical Sorption Models
Models which do not require state parameters. 


### Fundamental Sorption Models
Models which are described by an equation of state such as Sanchez Lacombe or PC-SAFT. 


### [Transient Sorption Models](ref)
Special models which describe time-dependent sorption into a flat polymer sheet (termed "slab") geometry. Useful for estimating diffusivity using measurements of sorption over time.


---

ainly accessed through the function `fit_model` and `predict_sorption`
