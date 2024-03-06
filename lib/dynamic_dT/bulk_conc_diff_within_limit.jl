function bulk_conc_diff_within_limit(new_bulk_concs, bulk_concs, constants_float)
    """
    This function determines whether all bulk concentrations changed less than the threshold 
    compared to the previous bulk concentrations

    Arguments
    new_bulk_concs:     A (ncompounds,) vector with the recently calculated bulk concentrations
    bulk_concs:         A (ncompounds,) vector with the previous bulk concentrations
    constants:          A "General" struct containing all the simulation constants

    Returns
    within_limit:       A Boolean indicating whether all bulk concentrations changed less than maximum permitted
    """

    within_limit = all(abs.(new_bulk_concs .- bulk_concs) ./ (bulk_concs .+ 1e-20) .< constants_float.maxRelDiffBulkConc)

    return within_limit
end
