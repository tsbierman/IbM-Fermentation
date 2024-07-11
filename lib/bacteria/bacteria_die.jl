function bacteria_die!(bac_vecfloat, bac_vecint, bac_vecbool, constants_float)
    """
    This function tests which bacteria are too small and sents them to the function killBacs!()

    Arguments
    bac_XYZ:            A struct containing bacterial parameters
    constants_float:    A "Float" struct containing simulation constants of type Float64

    Returns
    bac_XYZ:            A struct containing bacterial parameters without the bacteria that had to die
    """

    # Test whether the mass of bacteria is below minimum mass
    mask_tooSmall = bac_vecfloat.molarMass .* constants_float.bac_MW .< constants_float.min_bac_mass_grams
    nCellsTooSmall = sum(mask_tooSmall)

    # If any are too small --> go die
    if nCellsTooSmall != 0
        bac_vecfloat, bac_vecint, bac_vecbool = killBacs!(bac_vecfloat, bac_vecint, bac_vecbool, mask_tooSmall)
    end

    return bac_vecfloat, bac_vecint, bac_vecbool
end
