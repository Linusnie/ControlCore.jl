# discrete time conversions

convert{T}(::Type{DSisoZpk}, s::DSisoRational{T}) = zpk(s)
convert{T}(::Type{DSisoRational}, s::DSisoZpk{T}) = tf(s)

function tf{T}(s::DSisoZpk{T})
  DSisoRational(numpoly(s), denpoly(s), samplingtime(s))
end

function zpk{T}(s::DSisoRational{T})
  DSisoZpk(zpkdata(s)..., samplingtime(s))
end

# continuous time conversions

convert{T}(::Type{CSisoZpk}, s::CSisoRational{T}) = zpk(s)
convert{T}(::Type{CSisoRational}, s::CSisoZpk{T}) = tf(s)

function tf{T}(s::CSisoZpk{T})
  DSisoRational(numpoly(s), denpoly(s))
end

function zpk{T}(s::CSisoRational{T})
  DSisoZpk(zpkdata(s)...)
end
