function bacteria_shove!(bac_vecfloat, grid_float, grid_int, constants_float)
    """
    This function calls a Java function that shoves bacteria such that minimal overlap occurs.
    The Java function is used as it can use the QuadTree dstructure, which prevents that all
    bacteria are checked for distance and overlap.

    Arguments
    bac:                A "General" struct containing all parameters related to the bacteria
    grid:               A "General" struct containing all parameters related to the grid
    constants:          A "General" struct containing all the simulation constants

    Returns
    bac:                A bac struct with updates bacterial coordinates
    """

    # Import classes required
    BMQT_class = JavaCall.@jimport(shoving.BiomassQuadtree)
    Results_class = JavaCall.@jimport(shoving.Results)

    # Create instance of the QuadTree
    qt = BMQT_class((jdouble, jdouble, jdouble, jdouble), 0.0, grid_float.dx * grid_int.nx, 0.0, grid_float.dy * grid_int.ny)

    # Call shoving algorithm
    r = jcall(qt, "pushing2D", Results_class, (jint, Array{jdouble,1}, Array{jdouble,1}, Array{jdouble,1}, jdouble, jdouble, jdouble), length(bac_vecfloat.x), bac_vecfloat.x, bac_vecfloat.y, bac_vecfloat.radius, 0.1, constants_float.bac_max_radius * 2, constants_float.kDist)
    # Update the coordinates in struct
    bac_vecfloat.x = JavaCall.jfield(r, "bac_x", Array{jdouble,1})
    bac_vecfloat.y = JavaCall.jfield(r, "bac_y", Array{jdouble,1})
    return bac_vecfloat
end

# ----------------------- START VISUALISATION -----------------------------
# If some quick visualisation is desired, uncomment the following:
# using JavaCall
# BMQT_class = JavaCall.@jimport(shoving.BiomassQuadtree)
# Results_class = JavaCall.@jimport(shoving.Results)
# qt = BMQT_class((jdouble, jdouble, jdouble, jdouble), 0.0, grid_float.dx*grid_int.nx, 0.0, grid_float.dy*grid_int.ny)
# println(typeof(bac.x))
# r = jcall(qt, "pushing2D", Results_class, (jint, Array{jdouble,1}, Array{jdouble,1}, Array{jdouble,1}, jdouble, jdouble, jdouble), length(bac.x), bac.x, bac.y, bac.radius, 0.1, constants.bac_max_radius * 2, constants.kDist)
# bac_y = JavaCall.jfield(r, "bac_y", Array{jdouble,1})

# using Plots
# plot(bac.x, bac.y, seriestype=:scatter)
# ----------------------- START VISUALISATION -----------------------------
