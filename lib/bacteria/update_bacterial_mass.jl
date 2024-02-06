function update_bacterial_mass!(bac, dT)
    """
    This function integrates bacterial growth during the dT and updates the mass of each bacteria accordingly
    Depending on the growth rate (mu), the mass increases or decreases

    Arguments
    bac:                A "General" struct containing all parameters related to the bacteria
    dT:                 The timestep over which has to be integrated

    Returns
    bac:                A bac struct with updated bacterial mass
    """
    
    increase_index = bac.mu .>= 0 .& bac.active # Increase mass when mu positive and active bacteria
    decrease_index = bac.mu .<  0 .& bac.active # Decrease mass when mu negative and active bacteria
    
    # Updating
    bac.molarMass[increase_index] = bac.molarMass[increase_index] + dT * bac.mu[increase_index] .* bac.molarMass[increase_index]  
    bac.molarMass[decrease_index] = bac.molarMass[decrease_index] + dT * bac.mu[decrease_index] .* bac.molarMass[decrease_index]

    return bac
end
