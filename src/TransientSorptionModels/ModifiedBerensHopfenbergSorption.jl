struct ModifiedBerensHopfenbergSorption end

struct ModifiedBerensHopfenbergSorptionModel{MFT, KFT, MRT, KRT, BT} <: TransientSorptionModel
    m_f::MFT
    k_f::KFT
    m_r::MRT
    k_r::KRT
    beta::BT
    is_model_linearized::Bool
end

ModifiedBerensHopfenbergSorptionModel(mf, kf, mr, kr, beta) = ModifiedBerensHopfenbergSorptionModel(mf, kf, mr, kr, beta, false) 

function linearize_model(model::ModifiedBerensHopfenbergSorptionModel)
    if model.is_model_linearized
        throw(ErrorException("Model is already linearized"))
        return model
    else
        return ModifiedBerensHopfenbergSorptionModel(model.m_f, log(model.k_f), model.m_r, log(model.k_r), log(model.beta), true)
    end
end

function unlinearize_model(model::ModifiedBerensHopfenbergSorptionModel)
    if model.is_model_linearized
        return ModifiedBerensHopfenbergSorptionModel(model.m_f, exp(model.k_f), model.m_r, exp(model.k_r), exp(model.beta), false)
    else
        throw(ErrorException("Model is not linearized"))
        return model
    end
end

function predict_sorption(sorptionmodel::ModifiedBerensHopfenbergSorptionModel, time_seconds::Number; iters=10)
    if sorptionmodel.is_model_linearized
        sorptionmodel = unlinearize_model(sorptionmodel)
    end
    iter_start = 0
    if sorptionmodel.k_f/sorptionmodel.beta == 1
        iter_start = 1
    end  # if you look at the denominator, the zeroth term is undefined if these two are equal for the modified BH model
    fick_1 = 1

    beta_ratio = 4 * sorptionmodel.k_f / (pi^2 * sorptionmodel.beta)
    fick_2 = - exp(-sorptionmodel.beta * time_seconds) * sqrt(beta_ratio) * tan(sqrt(1/beta_ratio))

    # summation = sum([(exp(-(2*n+1)^2 * sorptionmodel.k_f * time_seconds)) / ((2*n+1)^2 * (1 - (2*n+1)^2 * sorptionmodel.k_f / sorptionmodel.beta)) for n in iter_start:(iters)])
    summation = 0 
    @inbounds for n in 0:iters
        summation += (exp(-(2*n+1)^2 * sorptionmodel.k_f * time_seconds)) / ((2*n+1)^2 * (1 - (2*n+1)^2 * sorptionmodel.k_f / sorptionmodel.beta))
    end
    
    fick_3 = -8/pi^2 * summation

    relax = sorptionmodel.m_r * (1 - exp(-sorptionmodel.k_r * time_seconds))
    
    total = sorptionmodel.m_f * (fick_1 + fick_2 + fick_3) + relax
    return total
end