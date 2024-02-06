function bacteria_divide!(bac, constants)
    """
    This function divides bacteria that are above the mass threshold. It splits the
    mass randomly between 45% and 55% and updates the radius.
    This is done in a while loop to make sure that even the already divided bacteria are not 
    above the threshold. If this is the case, then the timestep is too large.

    Arguments
    bac:                A "General" struct containing all parameters related to the bacteria
    constants:          A "General" struct containing all the simulation constants

    Returns
    bac:                A bac struct where every bacteria is below the maximum mass threshold.
    cycle:              The amount of cycles that had to be done to get every bacteria below the mass threshold
    """

    cycle = 0

    while sum(bac.molarMass * constants.bac_MW .> constants.max_bac_mass_grams) > 0
        cycle = cycle + 1

        # Select indices of too large bacteria
        mask_tooBig = bac.molarMass * constants.bac_MW .> constants.max_bac_mass_grams
        nCellsTooBig = sum(mask_tooBig)

        # Generate radial coordinates (random angle and fixed distance) and convert to cartesian coordiantes
        phi = rand(nCellsTooBig) * 2 * pi
        new_x = bac.x[mask_tooBig] + bac.radius[mask_tooBig] .* cos.(phi)
        new_y = bac.y[mask_tooBig] + bac.radius[mask_tooBig] .* sin.(phi)

        # Copy the other properties
        new_species = bac.species[mask_tooBig]
        new_mu = bac.mu[mask_tooBig]
        new_active = BitArray(ones(nCellsTooBig))

        # Split mass over parent and child
        new_molarMass = bac.molarMass[mask_tooBig] .* (0.45 .+ 0.1 * rand(nCellsTooBig))        # Mass of new cell is somewhere random between 0.45 and 0.55 of old mass
        bac.molarMass[mask_tooBig] = bac.molarMass[mask_tooBig] - new_molarMass                 # Parent cell keeps remaining part

        # Update radius of both
        new_radius = ((new_molarMass * constants.bac_MW / constants.bac_rho) * (3 / (4 * pi))) .^ (1/3)
        bac.radius[mask_tooBig] = ((bac.molarMass[mask_tooBig] * constants.bac_MW / constants.bac_rho) * (3 / (4 * pi))) .^ (1/3)

        # Update values
        bac.x = [bac.x; new_x]                                      # [-]
        bac.y = [bac.y; new_y]                                      # [-]
        bac.species = [bac.species; new_species]                    # [-]
        bac.molarMass = [bac.molarMass; new_molarMass]              # [mol]
        bac.radius = [bac.radius; new_radius]                       # [m]
        bac.mu = [bac.mu; new_mu]                                   # [h-1]
        bac.active = [bac.active; new_active]                       # [Boolean]
    end

    return bac, cycle
end
