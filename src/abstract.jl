abstract LtiSystem
abstract StateSpace       <: LtiSystem
abstract TransferFunction <: LtiSystem

abstract ContinuousTf     <: TransferFunction
abstract ContinuousSisoTf <: ContinuousTf

abstract DiscreteTf       <: TransferFunction
abstract DiscreteSisoTf   <: DiscreteTf

abstract ContinuousSs     <: StateSpace
abstract DiscreteSs       <: StateSpace

# TODO: When everything is finished, try to write catch-all functions with a
#       parametric style to see the difference, if any! Most likely, parametric
#       version is better (more efficient), i.e.,
#
#       +{T1<:LtiSystem, T2<:LtiSystem}(s1::T1, s2::T2) = +(promote(s1, s2)...)
+(sys1::LtiSystem, sys2::LtiSystem) = +(promote(sys1, sys2)...)
-(sys1::LtiSystem, sys2::LtiSystem) = -(promote(sys1, sys2)...)
*(sys1::LtiSystem, sys2::LtiSystem) = *(promote(sys1, sys2)...)
/(sys1::LtiSystem, sys2::LtiSystem) = /(promote(sys1, sys2)...)

series(sys1::LtiSystem, sys::LtiSystem)   = *(promote(sys1, sys2)...)
parallel(sys1::LtiSystem, sys::LtiSystem) = +(promote(sys1, sys2)...)
