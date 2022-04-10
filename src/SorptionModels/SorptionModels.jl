abstract type SorptionModel end


function predict_concentration end

"""
    fit_model(model, ::IsothermData)

Fit a sorption model to an isotherm. If the data required for the model is not in the isotherm, an error message will be returned. 

Any keyword arguments get passed on to the specific model.
"""
function fit_model(model, ::IsothermData) end