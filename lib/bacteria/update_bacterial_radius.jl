function update_bacterial_radius!(bac_vecfloat, constants)
    """
    This function updates the bacterial radius based on their masses.

    Arguments
    bac:                A "General" struct containing all parameters related to the bacteria
    constants:          A "General" struct containing all the simulation constants

    Returns
    bac:                A bac struct with updated radius
    """
    
    # Determine radius based on the mass and spherical cell
    bac_vecfloat.radius = ( (bac_vecfloat.molarMass * constants.bac_MW / constants.bac_rho) * (3 / (4 * pi)) ) .^ (1/3)
    return bac_vecfloat
end
