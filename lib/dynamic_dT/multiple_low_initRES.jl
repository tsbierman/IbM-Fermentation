function multiple_low_initRES(iProf, maxInitRES, Time, constants)
    """
    Determine whether the previous initial RES values were all below the threshold

    iProf:          Profiling Index (+1 for each simulated dT_bac)
    maxInitRES:     Maximum RES value of this RES iterations
    Time:           Struct containing all time variables
    constants:      Struct containing all system constants

    returns multiple_low: Boolean
    """

    recently_changed = (Time.current - Time.changed_dT) < constants.dynamicDT.nIterThreshold*Time.dT_bac # Checks if enough time has passed since last change of dT to increase dT
    if iProf >= constants.dynamicDT.nIterThreshold # can only check if iProf > n, otherwise negative indexing
        multiple_ss_with_low_init_RES = all(maxInitRES[iProf - constants.dynamicDT.nIterThreshold + 1 : iProf] .< constants.dynamicDT.initRESThresholdIncrease)
    else
        multiple_ss_with_low_init_RES = false
    end

    multiple_low = !recently_changed && multiple_ss_with_low_init_RES
    return multiple_low
end
