abstract LtiSystem
abstract StateSpace       <: LtiSystem
abstract TransferFunction <: LtiSystem

abstract ContinuousTf     <: TransferFunction
abstract DiscreteTf       <: TransferFunction

abstract ContinuousSisoTf <: ContinuousTf
abstract DiscreteSisoTf   <: DiscreteTf

+(sys1::LtiSystem, sys2::LtiSystem) = +(promote(sys1, sys2)...)
-(sys1::LtiSystem, sys2::LtiSystem) = -(promote(sys1, sys2)...)
*(sys1::LtiSystem, sys2::LtiSystem) = *(promote(sys1, sys2)...)
/(sys1::LtiSystem, sys2::LtiSystem) = /(promote(sys1, sys2)...)
