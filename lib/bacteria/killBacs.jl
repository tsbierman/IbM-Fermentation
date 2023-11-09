function killBacs!(bac, indices)
    """
    This function only selects and keeps the bacteria that are big enough using Inverse Indexing
    bac is a struct that contains all the information regarding the bacteria
    indices is a Vector of indices that indicate which bacteria are too small
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