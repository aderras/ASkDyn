#!/usr/bin/env julia1
push!(LOAD_PATH, string(pwd(),"/julia-files/"))

using UserInputs, Distributed, SpinDynamics

# Call the function that asks for user inputs
userParams = getuserinputs()

# Get confirmation
confirm = getconfirmation(userParams)
if confirm == false
    exit()
end

# If user selected parallel computing, launch cores
if userParams.llg.parallel == 1

  addprocs(userParams.llg.numCores)
  println("\nInitializing parallel cores. Number of processers = ",
    nprocs(),", number of workers = ", nworkers())

  # Need to define a module with the @everywhere macro so
  # every worker can run separately
  @everywhere module runScript

    # Set the load path for every worker to pwd()
    push!(LOAD_PATH, pwd())

    using UserInputs, SpinDynamics, Distributed

    function evaluate(p)

      allParams = []
      getparamlist!(p, allParams)

      pmap(launchcomputation, allParams)

    end

  end

  # Import the module into the current scope by using import
  # with a dot. Then call the parallel evaluation.
  using .runScript
  runScript.evaluate(userParams)

else

  # If the user did not request parallel evaluation, run
  # sequentially. Still need to get a list of all the
  # parameters in case user requested multiple values.
  allParams = []
  getparamlist!(userParams, allParams)

  map(launchcomputation, allParams)

end
