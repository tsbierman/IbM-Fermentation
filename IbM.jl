function IbM(simulation_number)
    """
    This is the main function to run the IbM model for the given preset file 
    (naming convention: sim_xxxx.mat, where xxxx is the simulation number)
    
    ->  creates a folder in the Results directory for the output of the simulation
    ->  Run the simulation based on the simulation file 
        (save results in the corresponding results folder)
    ->  Save profiling results in the corresponding results folder
    ->  Move the preset file in the corresponding results folder 
        (signifying that the simulation has been run) after running the model
    """

    # Get needed file
    include(string(pwd(), "\\lib\\integTime.jl"))

    # Argument check
    if !isa(simulation_number, Int) && !(0 <= simulation_number < 10000)
        throw(ArgumentError("Simulation number should be an integer in the range 1-9999"))
    end

    # Check if simulation file exists
    simulation_file = @sprintf("sim_%04d.jld2", simulation_number)
    if !isfile(simulation_file)
        throw(SystemError("The simulation file $(simulation_file) does not exist"))
    end
    
    # Create output directory for results
    output_dir = string(pwd(), @sprintf("\\results\\%04d", simulation_number))
    if !isdir(output_dir)
        mkdir(output_dir)
    end

    # Initialising the java should be done with initialiseJVM.jl beforehand
    # Making subfolders accesible is done with separate include statements (or play with push!(LOADPATH))

    # Dont know how to disable/enable specific warnings

    # ============== Time advancements ====================
    println("> SIMULATION RUNNING >>>>>>>>>> \n")
    tick()
    integTime(simulation_file, output_dir)
    totalTime = tok()

    println("> SIMULATION FINISHED >>>>>>>>>> \n")

    load(simulation_file, "constants")
    @sprrintf("\n\nTotal time for simulation of %.2f hours:\n\t%.2f seconds\n", constants.simulation_end, totalTime)

    # Cleanup of root directory
    mv(simulation_file, string(output_dir, "\\", simulation_file))

end
