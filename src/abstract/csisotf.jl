abstract CSisoTf{T} <: SisoTf{T}

isdiscrete(s::CSisoTf)          = false
samplingtime{T}(s::CSisoTf{T})  = zero(T)
