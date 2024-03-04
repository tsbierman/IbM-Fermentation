function multiple_high_iters(iDiffusion, iProf, nDiffIters, Time, constants_float, constants_vecint)
    """
    This function detects whether the previous steady states have been reached with
    high numbers of diffusion iterations.

    Arguments
    iDiffusion:         The number of cycles of diffusion (for this dT_bac)
    iProf:              The Profiling Index (+1 for each simulated dT_bac)
    nDiffIters:         A large storage for the number of diffusion iterations to reach steady state for each dT_bac
    Time:               A "General" struct containing all time-related parameters
    constants:          A "General" struct containing all the simulation constants

    Returns
    multiple_high:      A Boolean indicating whether the previous steady state have been reached with high numbers of diffusion iterations
    """

    changed_previous_ss = (Time.current - Time.changed_dT) < constants_vecint.nIterThreshold[1]*Time.dT_bac # Checks if enough time has passed since last change of dT to increase dT

    if iProf > constants_vecint.nIterThreshold[1]   # Check whether enough steadystates has passed since the start
        multiple_ss_with_high_nIter = all(vcat(nDiffIters[iProf-constants_vecint.nIterThreshold[1]+1:iProf-1], iDiffusion) .> constants_vecint.iterThresholdIncrease[1]) # Checks diffusion iterations with the threshold

    else
        multiple_ss_with_high_nIter = false         # Not enough time has passed, so false

    end

    multiple_high = !changed_previous_ss && multiple_ss_with_high_nIter

    return multiple_high
end
