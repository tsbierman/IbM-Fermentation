function bacteria_shove!(bac, grid, constants)
    """
    Functions calls a Java function that shoves bacteria so no overlap is present in the aggregate
    bac is a struct that contains all the information regarding the bacteria
    grid is a struct containing all infromation about the grid
    constants is a struct containing all simulation constants
    """

    BMQT_class = JavaCall.@jimport(shoving.BiomassQuadtree)
    Results_class = JavaCall.@jimport(shoving.Results)
    qt = BMQT_class((jdouble, jdouble, jdouble, jdouble), 0.0, grid.dx*grid.nx, 0.0, grid.dy*grid.ny)
    r = jcall(qt, "pushing2D", Results_class, (jint, Array{jdouble,1}, Array{jdouble,1}, Array{jdouble,1}, jdouble, jdouble, jdouble), length(bac.x), bac.x, bac.y, bac.radius, 0.1, constants.bac_max_radius * 2, constants.kDist)
    bac.x = JavaCall.jfield(r, "bac_x", Array{jdouble,1})
    bac.y = JavaCall.jfield(r, "bac_y", Array{jdouble,1})
    return bac
end

# If some quick visualisation is desired, uncomment the following:
# using JavaCall
# BMQT_class = JavaCall.@jimport(shoving.BiomassQuadtree)
# Results_class = JavaCall.@jimport(shoving.Results)
# qt = BMQT_class((jdouble, jdouble, jdouble, jdouble), 0.0, grid.dx*grid.nx, 0.0, grid.dy*grid.ny)
# println(typeof(bac.x))
# r = jcall(qt, "pushing2D", Results_class, (jint, Array{jdouble,1}, Array{jdouble,1}, Array{jdouble,1}, jdouble, jdouble, jdouble), length(bac.x), bac.x, bac.y, bac.radius, 0.1, constants.bac_max_radius * 2, constants.kDist)
# bac_y = JavaCall.jfield(r, "bac_y", Array{jdouble,1})

# using Plots
# plot(bac.x, bac.y, seriestype=:scatter)
