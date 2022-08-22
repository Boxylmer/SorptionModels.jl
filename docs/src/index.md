# SorptionModels.jl

SorptionModels.jl provides a number of models which describe solubility in polymers. They are broken down into three categories:

### Empirical Sorption Models
Models which do not require state parameters. 


### Fundamental Sorption Models
Models which are described by an equation of state such as Sanchez Lacombe or PC-SAFT. 


### [Transient Sorption Models](@ref)
Special models which describe time-dependent sorption into a flat polymer sheet (termed "slab") geometry. Useful for estimating diffusivity using measurements of sorption over time.


!!! warning
    Though all models are fit via `fit_model`, each particular model may require different types of inputs. 
    For example, the [Dual Mode](@ref) Model generally fits concentration to pressure or fugacity data, while the [GAB](@ref) Model fits concentration to activity. 

Ignoring the details of any specific model, all sorption model share the common aim to *predict penetrant-polymer solubility*. Therefore, every sorption model will have a few shared functions that can be documented here. Exceptions or features that exist but aren't implemented here should return a polite error message saying so.

Mainly accessed through the function `fit_model` and `predict_sorption`.
Generally, you can hand an isotherm (from MembraneBase) and the right supplementary information (if needed) over to `fit_model`.

```@docs
SorptionModels.fit_model(model, ::SorptionModels.IsothermData)
```

and then predict concentration via 
```@docs
SorptionModels.predict_concentration
```



