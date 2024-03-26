# Empirical Sorption Models

## Dual Mode

``C = k_D \cdot p + \frac{C_{H}^{'}bp}{1 + bp}``

The Dual Mode model is a semi-empirical model that combines Henry's Law (``k_D \cdot p``) with a *Langmuir-Hinshelwood*-like adsorption model (``\frac{C_{H}^{'}bp}{1 + bp}``) to describe sorption in glassy polymers. It has some predictive power and is widely used in the membrane field. 

```@autodocs
Modules = [SorptionModels]
Pages   = ["SorptionModels/DualMode.jl"]
```

## Flory-Huggins 

``C = φ / ((1-φ) \cdot V_{pen}``
where ``a = φ * exp((1-φ) + \chi * (1-φ)^2)``, ``a`` is activity and ``\chi`` is the Flory-Huggins interaction paramter. 

The Flory-Huggins model is a lattice fluid model that describes sorption via volume fraction and activity. This model is most typically used for rubbery polymers that exhibit swelling.

```@autodocs
Modules = [SorptionModels]
Pages   = ["SorptionModels/FloryHuggins.jl"]
```

## Flory-Huggins Dual Mode

``C = φ / ((1-φ) \cdot V_{pen} + \frac{C_{H}^{'}bp}{1 + bp}``
where ``a = φ * exp((1-φ) + \chi * (1-φ)^2)``, ``a`` is activity and ``\chi`` is the Flory-Huggins interaction paramter. 

This model combines the Dual Mode model with the Flory-Huggins model, essentially allowing for the description of glassy polymers which exhibit swelling.  
```@autodocs
Modules = [SorptionModels]
Pages   = ["SorptionModels/FloryHugginsDualMode.jl"]
```

## GAB

``C = \frac{C_pkAa}{(1-ka) \cdot (1 - ka - kAa)}``

The GAB model is a pure component model that describes polymer sorption as a function of activity.
It excels at vapor sorption, which typically demonstrates a sigmoidally shaped isotherm. Howerver, it is not very predictive, and is usually best suited for use as 
- an accurate interpolation method, or 
- to extract information from the isotherm in the form of the model paramters, which have physical meaning. 

```@autodocs
Modules = [SorptionModels]
Pages   = ["SorptionModels/GAB.jl"]
```