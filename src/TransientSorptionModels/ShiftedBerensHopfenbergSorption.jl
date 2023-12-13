struct ShiftedBerensHopfenbergSorption end

struct ShiftedBerensHopfenbergSorptionModel{MFT, KFT, MRT, KRT, TST} <: TransientSorptionModel
    m_f::MFT
    k_f::KFT
    m_r::MRT
    k_r::KRT
    t_s::TST
    is_model_linearized::Bool
end

ShiftedBerensHopfenbergSorptionModel(mf, kf, mr, kr, ts) = ShiftedBerensHopfenbergSorptionModel(mf, kf, mr, kr, ts, false) 

function linearize_model(model::ShiftedBerensHopfenbergSorptionModel)
    if model.is_model_linearized
        throw(ErrorException("Model is already linearized."))
        return model
    else
        return ShiftedBerensHopfenbergSorptionModel(model.m_f, log(model.k_f), model.m_r, log(model.k_r), model.t_s, true)
    end
end

function unlinearize_model(model::ShiftedBerensHopfenbergSorptionModel)
    if model.is_model_linearized
        return ShiftedBerensHopfenbergSorptionModel(model.m_f, exp(model.k_f), model.m_r, exp(model.k_r), model.t_s, false)
    else
        throw(ErrorException("Model was not linearized."))
        return model
    end
end

function predict_sorption(sorptionmodel::ShiftedBerensHopfenbergSorptionModel, time_seconds::Number; iters=15)
    if sorptionmodel.is_model_linearized
        _sorptionmodel = unlinearize_model(sorptionmodel)
    else
        _sorptionmodel = sorptionmodel
    end
    summation = 0
    @inbounds for n in 0:iters
        summation += 1/(2*n+1)^2 * exp(-(2*n+1)^2 * _sorptionmodel.k_f * (time_seconds + model.t_s))
    end

    _sorptionmodel.m_f * (1 - 8/(pi^2) * summation) + _sorptionmodel.m_r * (1 - exp(-_sorptionmodel.k_r * (time_seconds + model.t_s)))
end

