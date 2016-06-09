abstract LtiSystem

# Printing functions
summary(s::LtiSystem) = string("lti(nx=", numstates(s), ",ny=", numoutputs(s),
  ",nu=", numinputs(s), (isdiscrete(s) ? ",Ts=", samplingtime(s), : ""), ")")

showcompact(io::IO, s::LtiSystem) = print(io, summary(s))
show(io::IO, s::LtiSystem)        = print(io, summary(s))
showall(io::IO, s::LtiSystem)     = print(io, summary(s))
