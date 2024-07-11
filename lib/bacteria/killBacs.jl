function killBacs!(bac_vecfloat, bac_vecint, bac_vecbool, indices)
    """
    This function removes the indices corresponding to organisms that did not fullfill requirements.
    The indices that did fullfill them are selected with Inverse Indexing

    Arguments
    bac_XYZ:            A struct containing bacterial parameters
    indices:            A vector which indicates which bacteria have to be killed

    Returns
    bac_XYZ:            A struct containing bacterial parameters, but without the bacterial that had to be killed
    """
    bac_vecfloat.x = bac_vecfloat.x[Not(indices)]
    bac_vecfloat.y = bac_vecfloat.y[Not(indices)]
    bac_vecfloat.radius = bac_vecfloat.radius[Not(indices)]
    bac_vecint.species = bac_vecint.species[Not(indices)]
    bac_vecfloat.molarMass = bac_vecfloat.molarMass[Not(indices)]
    bac_vecbool.active = bac_vecbool.active[Not(indices)]
    bac_vecfloat.mu = bac_vecfloat.mu[Not(indices)]
    bac_vecint.colony_nums = bac_vecint.colony_nums[Not(indices)]
    return bac_vecfloat, bac_vecint, bac_vecbool
end
