abstract DSisoTf{T} <: SisoTf{T}

isdiscrete(s::DSisoTf)          = true
samplingtime{T}(s::DSisoTf{T})  = s.Ts
