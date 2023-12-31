function granule_density!(bac, constants)
    """
    Updates the density of the granule
    bac is a struct containing all bacteria information
    constants is a struct containing simulation constants
    """
    bac_m = bac.molarMass * constants.bac_MW                       # [g] for every bacteria
    # Approximated as a square due to irregular growth directions
    bac.bac_rho_bio = sum(bac_m) / ((maximum(bac.y) - minimum(bac.y)) * (maximum(bac.x) - minimum(bac.x)) * 1e-6) # [g/m3] The 1e-6 is the depth dimension
    println(bac.bac_rho_bio)
    return bac
end