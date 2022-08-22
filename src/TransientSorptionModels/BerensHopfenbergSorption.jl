struct BerensHopfenbergSorption end

struct BerensHopfenbergSorptionModel{MFT, KFT, MRT, KRT} <: TransientSorptionModel
    m_f::MFT
    k_f::KFT
    m_r::MRT
    k_r::KRT
    is_model_linearized::Bool
end

BerensHopfenbergSorptionModel(mf, kf, mr, kr) = BerensHopfenbergSorptionModel(mf, kf, mr, kr, false) 

function linearize_model(model::BerensHopfenbergSorptionModel)
    if model.is_model_linearized
        throw(ErrorException("Nope"))
        return model
    else
        return BerensHopfenbergSorptionModel(model.m_f, log(model.k_f), model.m_r, log(model.k_r), true)
    end
end

function unlinearize_model(model::BerensHopfenbergSorptionModel)
    if model.is_model_linearized
        return BerensHopfenbergSorptionModel(model.m_f, exp(model.k_f), model.m_r, exp(model.k_r), false)
    else
        throw(ErrorException("Nope"))
        return model
    end
end

function predict_sorption(sorptionmodel::BerensHopfenbergSorptionModel, time_seconds::Number; iters=15)
    if sorptionmodel.is_model_linearized
        _sorptionmodel = unlinearize_model(sorptionmodel)
    else
        _sorptionmodel = sorptionmodel
    end
    # summation = sum([1/(2*n+1)^2 * exp(-(2*n+1)^2 * sorptionmodel.k_f * time_seconds) for n in 0:(iters)])
    summation = 0
    @inbounds for n in 0:iters
        summation += 1/(2*n+1)^2 * exp(-(2*n+1)^2 * _sorptionmodel.k_f * time_seconds)
    end

    _sorptionmodel.m_f * (1 - 8/(pi^2) * summation) + _sorptionmodel.m_r * (1 - exp(-_sorptionmodel.k_r * time_seconds))
end

