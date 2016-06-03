abstract LtiSystem
abstract StateSpace       <: LtiSystem
abstract TransferFunction <: LtiSystem

abstract ContinuousTf     <: TransferFunction
abstract ContinuousSisoTf <: ContinuousTf

abstract DiscreteTf       <: TransferFunction
abstract DiscreteSisoTf   <: DiscreteTf

abstract ContinuousSs     <: StateSpace
abstract DiscreteSs       <: StateSpace

+(sys1::LtiSystem, sys2::LtiSystem) = +(promote(sys1, sys2)...)
-(sys1::LtiSystem, sys2::LtiSystem) = -(promote(sys1, sys2)...)
*(sys1::LtiSystem, sys2::LtiSystem) = *(promote(sys1, sys2)...)
/(sys1::LtiSystem, sys2::LtiSystem) = /(promote(sys1, sys2)...)
