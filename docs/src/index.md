# SorptionModels.jl

*Get more out of your solubility data*

SorptionModels.jl provides a number of models which describe solubility in polymers, as well as various analyses you can perform to extract information from your models once they are fit to experimental data.


## Using this package
Most methods are designed to operate on **isotherms**. This is accessed, computationally speaking, through `IsothermData` structs, defined in MemrbaneBase. # todo link this here

!!! tip 
    Some models have convenience methods that circumvent the need to wrap data in an isotherm, but you are more likely to encounter user mistakes by doing so. For example, when fitting [Dual Mode](@ref) models for use in a [Partial Immobilization](@ref) analysis, you must be careful to ensure you're using fugacity instead of pressure. Ensuring that you're fitting with fugacities is as easy as setting the keyword `use_fugacity=true` when using IsothermStructs.

Models are broken down into three categories, each serving a particular range of purposes with some overlap:

### [Empirical Sorption Models](@ref)
Models which do not require state parameters. The parameters you fit from these models will typically have physical meaning which can be useful in understanding how sorption behaves in the polymer, but the models themselves are generally not incredibly predictive unless heavily leveraged against carefully taken experimental data (such as multiple isotherms at distinct temperatures). These models can also serve as accurate interpolation (and possibly extrapolation / noise suppression) tools. 


### [Fundamental Sorption Models](@ref)
Models which are described by an equation of state such as Sanchez Lacombe or PC-SAFT. These models, in stark contrast to empirical models, are predictive by nature. That is to say, they will generally not align *perfectly* to experimental data, but will generate values without much (if any) experimental data to go off of. They are also useful for predicting the solubility of mixed systems. 

To learn how to set up the equations of state that act as the engine of these models, see SorptionModels.jl.


### [Transient Sorption Models](@ref)
Special models which describe time-dependent sorption into a flat polymer sheet (termed "slab") geometry. Useful for estimating diffusivity using measurements of sorption over time (termed "sorption kinetics"), if such data is available. 

---

!!! warning
    Though all models are fit via `fit_model`, each particular model may require different types of inputs. 
    For example, the [Dual Mode](@ref) Model generally fits concentration to pressure or fugacity data, while the [GAB](@ref) Model fits concentration to activity. 

Obviously, all sorption models share the common aim to *predict penetrant-polymer solubility*. Therefore, every sorption model will have a few shared functions that can be understood without knowing individual model details.

All models are fit via `fit_model` and generate data through `predict_sorption`.
Generally, you can hand an isotherm (from MembraneBase) and the right supplementary information (if needed) over to... 
```@docs
SorptionModels.fit_model(model, ::SorptionModels.IsothermData)
```
...then immediately get the fitted isotherm back through...
```@docs
SorptionModels.predict_concentration
```

---



