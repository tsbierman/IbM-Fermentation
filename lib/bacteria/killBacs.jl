function killBacs!(bac, indices)
    """
    This function selects the indices it needs to keep by using Inverse Indexing

    Arguments
    bac:                A "General" struct containing all parameters related to the bacteria
    indices:            The indices of the bacteria that have to be killed

    Returns
    bac                 A bac struct with the bacteria that had to be killed removed
    """
    bac.x = bac.x[Not(indices)]
    bac.y = bac.y[Not(indices)]
    bac.radius = bac.radius[Not(indices)]
    bac.species = bac.species[Not(indices)]
    bac.molarMass = bac.molarMass[Not(indices)]
    bac.active = bac.active[Not(indices)]
    bac.mu = bac.mu[Not(indices)]
    return bac
end
