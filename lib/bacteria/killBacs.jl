function killBacs!(bac_vecfloat, bac_vecint, bac_vecbool, indices)
    """
    This function selects the indices it needs to keep by using Inverse Indexing

    Arguments
    bac:                A "General" struct containing all parameters related to the bacteria
    indices:            The indices of the bacteria that have to be killed

    Returns
    bac                 A bac struct with the bacteria that had to be killed removed
    """
    bac_vecfloat.x = bac_vecfloat.x[Not(indices)]
    bac_vecfloat.y = bac_vecfloat.y[Not(indices)]
    bac_vecfloat.radius = bac_vecfloat.radius[Not(indices)]
    bac_vecint.species = bac_vecint.species[Not(indices)]
    bac_vecfloat.molarMass = bac_vecfloat.molarMass[Not(indices)]
    bac_vecbool.active = bac_vecbool.active[Not(indices)]
    bac_vecfloat.mu = bac_vecfloat.mu[Not(indices)]
    return bac_vecfloat, bac_vecint, bac_vecbool
end
