# SorptionModels

[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://Boxylmer.github.io/SorptionModels.jl/dev/)
[![Build Status](https://github.com/Boxylmer/SorptionModels.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/Boxylmer/SorptionModels.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Coverage](https://codecov.io/gh/Boxylmer/SorptionModels.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/Boxylmer/SorptionModels.jl/tree/master)


SorptionModels is a package that helps you make the most of your solubility data. 

Quickstart: 

```
using SorptionModels
using MembraneBase

pressures = [0, 0.054409635, 0.101117776, 0.17570818, 0.275518841, 0.373215869]
concentrations = [0, 5.084334416, 8.57780684, 12.92845101, 17.55678063, 21.17207162]
my_isotherm = IsothermData(
    partial_pressures_mpa = pressures, 
    concentrations_cc = concentrations
)

dualmode_model = fit_model(DualMode(), my_isotherm)
predicted = predict_concentration(dualmode_model, pressures)

@show dualmode_model
@show concentrations
@show predicted
```