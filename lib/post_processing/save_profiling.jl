function save_profiling(profiling, maxErrors, normOverTime, nDiffIters, bulk_history, Time, directory)
    """
    This function saves important profiling arrays with information on the performance
    of the simulation. It will overwrite the last backup in order to always 
    have the latest file.
    
    profiling:          A matrix with per dT_bac the time per functionality of the model
    maxErrors:          A vector with per dT_bac the maximum RES value
    normOverTime:       A vector with per dT_bac the norm of delta-concentrations
    nDiffIters:         A vector with per dT_bac the number of diffusion iterations
    bulk_history:       A vector with per dT_bac the bulk concentration of each compound
    Time                A Time struct with decreased dT (diffusion)
    directory:          The directory where results are to be stored in
    """

    # Name file
    results_file = string(directory, "\\profilingResults.jdl2")

    # Save everything
    save(results_file, "profiling", profiling, "maxErrors", maxErrors, "normOverTime", normOverTime, "nDiffIters", nDiffIters, "bulk_history", bulk_history, "Time", Time)
end
