function bulk_conc_diff_within_limit(new_bulk_concs, bulk_concs, constants)
    """
    Determine whether all bulk concentrations changed less than the threshold compared to
    the previous bulk concentrations
    """

    within_limit = all(abs.(new_bulk_concs .- bulk_concs) ./ (bulk_concs .+ 1e-20) .< constants.dynamicDT.maxRelDiffBulkConc)

    return within_limit
end
