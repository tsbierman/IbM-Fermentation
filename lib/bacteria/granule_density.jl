function granule_density!(bac_vecfloat, constants)
    """
    Updates the density of the granule

    Arguments
    bac:                A "General" struct containing all parameters related to the bacteria
    constants:          A "General" struct containing all the simulation constants

    Returns
    bac:                A bac struct with an updated biological density of the granule
    """

    bac_m = bac_vecfloat.molarMass * constants.bac_MW                                                                        # [g] for every bacteria

    # Approximated as a square due to irregular growth directions
    bac_vecfloat.bac_rho_bio = sum(bac_m) / ((maximum(bac_vecfloat.y) - minimum(bac_vecfloat.y)) * (maximum(bac_vecfloat.x) - minimum(bac_vecfloat.x)) * 2e-6)   # [g/m3] The 2e-6 is the depth dimension
    return bac_vecfloat
end
