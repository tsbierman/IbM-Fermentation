function bacteria_shove(bac, grid, constants)
    """
    Functions calls a Java function that shoves bacteria so no overlap is present in the aggregate
    bac is a struct that contains all the information regarding the bacteria
    constants is a struct containing all simulation constants
    """
    JavaCall.init()
    jlm = @jimport "BiomassQuadtree"

    qt = shoving.BiomassQuadtree(0, grid.dx * grid.nx, 0, grid.dy*grid.ny)
    r = qt.pushing2D(length(bac.x), bac.x, bac.y, bac.radius, 0.1, constants.bac_max_radius * 2, constants.kDist)
    bac.x = r.bac_x
    bac.y = r.bac_y
    return bac
end
bacteria_shove(1,1,1)