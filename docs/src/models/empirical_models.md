# Empirical Sorption Models

## Dual Mode

``C = k_D \cdot p + \frac{C_{H}^{'}bp}{1 + bp}``

The Dual Mode model is a semi-empirical model that combines Henry's Law (``k_D \cdot p``) with a *Langmuir-Hinshelwood*-like adsorption model (``\frac{C_{H}^{'}bp}{1 + bp}``) to describe sorption in glassy polymers. It has some predictive power and is widely used in the membrane field. 

```@docs
SorptionModels.fit_model(::DualMode, ::SorptionModels.IsothermData; uncertainty_method)
SorptionModels.predict_concentration(::DualModeModel, pressures_mpa::AbstractVector)
SorptionModels.predict_concentration(::DualModeModel, pressure_mpa::Number)
SorptionModels.predict_pressure(::DualModeModel, ::Number)
```


## GAB

``C = \frac{C_pkAa}{(1-ka) \cdot (1 - ka - kAa)}``

The GAB model is a pure component model that describes polymer sorption as a function of activity.
It excels at vapor sorption, which typically demonstrates a sigmoidally shaped isotherm. Howerver, it is not very predictive, and is usually best suited for use as 
- an accurate interpolation method, or 
- to extract information from the isotherm in the form of the model paramters, which have physical meaning. 

```@docs
SorptionModels.fit_model(::GABModel, isotherm::SorptionModels.IsothermData; kwargs...)
SorptionModels.fit_gab_model(activities::AbstractVector, concentrations::AbstractVector; uncertainty_method=nothing)
SorptionModels.predict_concentration(gm::GABModel, activities::AbstractVector)
SorptionModels.predict_concentration(gm::GABModel, activity::Number)
```
