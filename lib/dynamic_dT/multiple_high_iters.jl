function multiple_high_iters(iDiffusion, iProf, nDiffIters, Time, constants)
    """
    Detect whether the previous steady states have been reached with
    a high number of diffusion iterations.
    iDiffusion:     number of cycles of diffusion (for this dT_bac)
    iProf:          Profiling Index (+1 for each simulated dT_bac)
    nDiffIters:     Storage for number of diffusion iterations to reach steady state for each dT_bac
    Time:           Struct containing all time variables
    constants:      Struct containing all system constants
    """

    changed_previous_ss = (Time.current - Time.changed_dT) < constants.dynamicDT.nIterThreshold*Time.dT_bac # Checks if enough time has passed since last change of dT to increase dT
    if iProf > constants.dynamicDT.nIterThreshold # Check whether enough steadystates has passed since the start
        multiple_ss_with_high_nIter = all(vcat(nDiffIters[iProf-constants.dynamicDT.nIterThreshold+1:iProf-1], iDiffusion) .> constants.dynamicDT.iterThresholdIncrease)
    else
        multiple_ss_with_high_nIter = false
    end

    multiple_high = !changed_previous_ss && multiple_ss_with_high_nIter

    return multiple_high
end
