#!/usr/bin/env julia1
push!(LOAD_PATH, pwd())

module main_run

  using initialCondition, userInputs, spinDynamics

  params = getUserParams()
  # params = getTestParams()

  # # println( params.mp )
  # println( params.llg )
  # println( params.fa )
  # println( params.ic )
  # println( params.pin )
  # println( params.defect )
  # println( params.save )

  s0 = buildInitial(params.ic, params.mp)

  # println( [getfield( params.mp,x ) for x in fieldnames( typeof(params.mp) ) ] )
  # println( calcEnergy(s0, params.mp, params.defect) )
  evaluateLL!( s0, params )

end
