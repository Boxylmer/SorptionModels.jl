struct FickianSorptionModel{MFT, KFT} <: TransientSorptionModel
    m_f::MFT
    k_f::KFT
    is_model_linearized::Bool
end

FickianSorptionModel(mf, kf) = FickianSorptionModel(mf, kf, false) 

function linearize_model(model::FickianSorptionModel)
    if model.is_model_linearized
        throw(ErrorException("Nope"))
        return model
    else
        return FickianSorptionModel(model.m_f, log(model.k_f), true)
    end
end

function unlinearize_model(model::FickianSorptionModel)
    if model.is_model_linearized
        return FickianSorptionModel(model.m_f, exp(model.k_f), false)
    else
        throw(ErrorException("Nope"))
        return model
    end
end

function predict_sorption(sorptionmodel::FickianSorptionModel, time_seconds::Number; iters=15)
    if sorptionmodel.is_model_linearized
        sorptionmodel = unlinearize_model(sorptionmodel)
    end
    summation = sum([1/(2*n+1)^2 * exp(-(2*n+1)^2 * sorptionmodel.k_f * time_seconds) for n in 0:(iters)])
    sorptionmodel.m_f * (1 - 8/pi^2 * summation)
end

