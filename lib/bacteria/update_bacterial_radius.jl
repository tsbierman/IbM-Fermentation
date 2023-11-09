function update_bacterial_radius!(bac, constants)
    """
    This function updates the radius of all bacteria
    bac is a struct containing all information regarding the bacteria
    constants is a struct containing all simulation constants
    """
    # Determine radius based on the mass and spherical cell
    bac.radius = ( (bac.molarMass * constants.bac_MW / constants.bac_rho) * (3 / (4 * pi)) ) .^ (1/3)
    return bac
end
