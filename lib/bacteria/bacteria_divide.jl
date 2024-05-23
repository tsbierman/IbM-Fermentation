function bacteria_divide!(bac_vecfloat, bac_vecint, bac_vecbool, constants_float)
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

    while sum(bac_vecfloat.molarMass .* constants_float.bac_MW .> constants_float.max_bac_mass_grams) > 0
        cycle = cycle + 1

        # Select indices of too large bacteria
        mask_tooBig = bac_vecfloat.molarMass .* constants_float.bac_MW .> constants_float.max_bac_mass_grams
        nCellsTooBig = sum(mask_tooBig)

        # Generate radial coordinates (random angle and fixed distance) and convert to cartesian coordiantes
        phi = rand(nCellsTooBig) .* 2 .* pi
        new_x = bac_vecfloat.x[mask_tooBig] .+ bac_vecfloat.radius[mask_tooBig] .* cos.(phi)
        new_y = bac_vecfloat.y[mask_tooBig] .+ bac_vecfloat.radius[mask_tooBig] .* sin.(phi)

        # Copy the other properties
        new_species = bac_vecint.species[mask_tooBig]
        new_mu = bac_vecfloat.mu[mask_tooBig]
        new_active = BitArray(ones(nCellsTooBig))
        new_colony_nums = bac_vecint.colony_nums[mask_tooBig]

        # Split mass over parent and child
        new_molarMass = bac_vecfloat.molarMass[mask_tooBig] .* (0.45 .+ 0.1 .* rand(nCellsTooBig))        # Mass of new cell is somewhere random between 0.45 and 0.55 of old mass
        bac_vecfloat.molarMass[mask_tooBig] = bac_vecfloat.molarMass[mask_tooBig] .- new_molarMass                 # Parent cell keeps remaining part

        # Update radius of both
        new_radius = ((new_molarMass .* constants_float.bac_MW ./ constants_float.bac_rho) .* (3 ./ (4 .* pi))) .^ (1 ./ 3)
        bac_vecfloat.radius[mask_tooBig] = ((bac_vecfloat.molarMass[mask_tooBig] .* constants_float.bac_MW ./ constants_float.bac_rho) .* (3 ./ (4 .* pi))) .^ (1 ./ 3)

        # Update values
        bac_vecfloat.x = [bac_vecfloat.x; new_x]                                      # [-]
        bac_vecfloat.y = [bac_vecfloat.y; new_y]                                      # [-]
        bac_vecint.species = [bac_vecint.species; new_species]                    # [-]
        bac_vecfloat.molarMass = [bac_vecfloat.molarMass; new_molarMass]              # [mol]
        bac_vecfloat.radius = [bac_vecfloat.radius; new_radius]                       # [m]
        bac_vecfloat.mu = [bac_vecfloat.mu; new_mu]                                   # [h-1]
        bac_vecbool.active = [bac_vecbool.active; new_active]                       # [Boolean]
        bac_vecint.colony_nums = [bac_vecint.colony_nums; new_colony_nums]
    end

    return bac_vecfloat, bac_vecint, bac_vecbool, cycle
end
