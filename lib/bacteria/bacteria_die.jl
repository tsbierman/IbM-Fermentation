function bacteria_die!(bac, constants)
    """
    This function tests which bacteria are too small and sents them to die
    bac is a struct that contains all the information regarding the bacteria
    constants is a struct containing all simulation constants
    """
    include("killBacs.jl")

    mask_tooSmall = bac.molarMass * constants.bac_MW .< constants.min_bac_mass_grams
    nCellsTooSmall = sum(mask_tooSmall)

    if nCellsTooSmall != 0
        bac = killBacs!(bac, mask_tooSmall)
    end

    return bac
end
