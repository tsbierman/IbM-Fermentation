function multiple_low_initRES(iProf, maxInitRES, Time, constants)
    """
    This function determines whether the previous initial RES values were all below the threshold

    Arguments
    iProf:              A Profiling Index (+1 for each simulated dT_bac)
    maxInitRES:         The maximum RES value of this RES iterations
    Time:               A "General" struct containing all time-related parameters
    constants:          A "General" struct containing all the simulation constants

    Returns 
    multiple_low:       Boolean indicating whether the previous initial RES values were below threshold
    """

    recently_changed = (Time.current - Time.changed_dT_bac) < constants.dynamicDT.nIterThreshold*Time.dT_bac # Checks if enough time has passed since last change of dT to increase dT

    if iProf >= constants.dynamicDT.nIterThreshold  # Can only check if iProf > n, otherwise negative indexing
        multiple_ss_with_low_init_RES = all(maxInitRES[iProf - constants.dynamicDT.nIterThreshold + 1 : iProf] .< constants.dynamicDT.initRESThresholdIncrease) # Checks RES against threshold

    else
        multiple_ss_with_low_init_RES = false       # Not enough time has passed, so false

    end

    multiple_low = !recently_changed && multiple_ss_with_low_init_RES
    
    return multiple_low
end
