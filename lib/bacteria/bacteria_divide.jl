function bacteria_divide!(bac, constants)
    """
    Divide bacteria that are above mass threshold
    bac is a struct containing all information regarding the bacteria
    constants is a struct containing all simulation constants

    The while loop is too make sure that the bacteria are not so big they need
    to divide twice in a row.
    """

    cycle = 0

    while sum(bac.molarMass * constants.bac_MW .> constants.max_bac_mass_grams) > 0
        cycle = cycle + 1

        # Select indices of too large bacteria
        mask_tooBig = bac.molarMass * constants.bac_MW .> constants.max_bac_mass_grams
        nCellsTooBig = sum(mask_tooBig)

        # Generate radial coordinates (random angle, fixed distance) and convert to cartesian coordiantes
        phi = rand(nCellsTooBig) * 2 * pi
        new_x = bac.x[mask_tooBig] + bac.radius[mask_tooBig] .* cos(phi)
        new_y = bac.y[mask_tooBig] + bac.radius[mask_tooBig] .* sin(phi)

        # Copy the other properties
        new_species = bac.species[mask_tooBig]
        new_mu = bac.mu[mask_tooBig]
        new_active = BitArray(ones(nCellsTooBig))

        # Split mass over parent and child
        new_molarMass = bac.molarMass[mask_tooBig] .* (0.45 + 0.1 * rand(nCellsTooBig)) # Mass of new cell is somewhere random between 0.45 and 0.55 of old mass
        bac.molarMass[mask_tooBig] = bac.molarMass[mask_tooBig] - new_molarMass # Parent cell keeps remaining part
        # Update radius of both
        new_radius = ((new_molarMass * constants.bac_MW / constants.bac_rho) * (3 / (4 * pi))) .^ (1/3)
        bac.radius[mask_tooBig] = ((bac.molarMass[mask_tooBig] * constants.bac_MW / constants.bac_rho) * (3 / (4 * pi))) .^ (1/3)

        # Update values
        bac.x = [bac.x; new_x]
        bac.y = [bac.y; new_y]
        bac.species = [bac.species; new_species]
        bac.molarMass = [bac.molarMass; new_molarMass]
        bac.radius = [bac.radius; new_radius]
        bac.mu = [bac.mu; new_mu]
        bac.active = [bac.active; new_active]
    end

    return bac, cycle
end
