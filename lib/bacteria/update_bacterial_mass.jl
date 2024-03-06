function update_bacterial_mass!(bac_vecfloat, bac_vecbool, dT)
    """
    This function integrates bacterial growth during the dT and updates the mass of each bacteria accordingly
    Depending on the growth rate (mu), the mass increases or decreases

    Arguments
    bac:                A "General" struct containing all parameters related to the bacteria
    dT:                 The timestep over which has to be integrated

    Returns
    bac:                A bac struct with updated bacterial mass
    """
    
    increase_index = bac_vecfloat.mu .>= 0 .& bac_vecbool.active # Increase mass when mu positive and active bacteria
    decrease_index = bac_vecfloat.mu .<  0 .& bac_vecbool.active # Decrease mass when mu negative and active bacteria
    
    # Updating
    bac_vecfloat.molarMass[increase_index] = bac_vecfloat.molarMass[increase_index] + dT * bac_vecfloat.mu[increase_index] .* bac_vecfloat.molarMass[increase_index]  
    bac_vecfloat.molarMass[decrease_index] = bac_vecfloat.molarMass[decrease_index] + dT * bac_vecfloat.mu[decrease_index] .* bac_vecfloat.molarMass[decrease_index]

    return bac_vecfloat
end
