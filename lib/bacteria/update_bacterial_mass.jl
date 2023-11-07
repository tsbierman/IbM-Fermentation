function update_bacterial_mass!(bac, dT)
    """
    This function integrates the bacterial growth during the dT and updates the mass accordingly
    bac is a struct containing all information regarding the bacteria
    dT is the time over which to integrates
    """
    
    increase_index = bac.mu .>= 0 .& bac.active # Increase mass when mu positive and active bacteria
    decrease_index = bac.mu .<  0 .& bac.active # Decrease mass when mu negative and active bacteria
    
    # Updating
    bac.molarMass[increase_index] = bac.molarMass[increase_index] + dT * bac.mu[increase_index] .* bac.molarMass[increase_index]  
    bac.molarMass[decrease_index] = bac.molarMass[decrease_index] + dT * bac.mu[decrease_index] .* bac.molarMass[decrease_index]

    return bac
end
