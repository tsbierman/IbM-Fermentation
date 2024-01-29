function bacteria_die!(bac, constants)
    """
    This function tests which bacteria are too small and sents them to the function killBacs!()

    Arguments
    bac:                A "General" struct containing all parameters related to the bacteria
    constants:          A "General" struct containing all the simulation constants

    Returns
    bac:                A bac struct without the bacteria that had to die
    """
    # Include required file
    # include(string(pwd(), "\\lib\\bacteria\\killBacs.jl"))

    # Test whether the mass of bacteria is below minimum mass
    mask_tooSmall = bac.molarMass * constants.bac_MW .< constants.min_bac_mass_grams
    nCellsTooSmall = sum(mask_tooSmall)

    # If any are too small --> go die
    if nCellsTooSmall != 0
        bac = killBacs!(bac, mask_tooSmall)
    end

    return bac
end
