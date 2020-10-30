#!/usr/bin/env julia1
push!(LOAD_PATH, pwd())

using userInputs, Distributed, spinDynamics

testParams = getTestParams()

println( "Test params mp: ", testParams.mp )
println( "Test params llg: ", testParams.llg )
println( "Param ranges: ", testParams.range )

if testParams.llg.parallel == 1

  addprocs( testParams.llg.numCores )
  println( "Initializing parallel cores. Number of processers = ", nprocs(),", number of workers = ", nworkers() )

    @everywhere module Test
    push!(LOAD_PATH, pwd())

    using userInputs, spinDynamics, Distributed

    function dosomething(testParams)

      allParams = []
      getParamList!( testParams, allParams )

      pmap( evaluateSpinDynamics, allParams )

    end

  end

  using .Test 
  Test.dosomething(testParams)

else 

  allParams = []
  getParamList!( testParams, allParams )

  map( evaluateSpinDynamics, allParams )

end



