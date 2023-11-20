function calcTimeOfDetach(bac, grid, grid2bac, grid2nBacs, constants)
    """
    This function calculates the time of detachment for each gridcell in the simulation domain
    Bac is a struct that contains all information regarding the bacteria
    grid is a struct containing all information regarding spacial discretization
    grid2bac is a nx * nY * ? matrix that for each cell contains which bacteria is located 
    there. The number corresponds to the index in the bac struct
    grid2nBacs is a (nx * ny) matrix containing the number of bacteria in each grid cell

    It returns a matrix with per gridcell the time of detachment
    """
    # Extract variables
    kDet = constants.kDet

    # Calculate centre of granule
    x_centre = mean(bac.x)
    y_centre = mean(bac.y)

    # Where are bacteria located in the grid? The grid2nBacs is slightly extended and 
    # morphed into a logical matrix
    detachment_grid = deepcopy(grid) # To keep the original grid values from changing
    detachment_grid.blayer_thickness = constants.kDist * constants.bac_max_radius * 2 # Distance factor times diameter

    # Returns logical matrix with per cell whether it is in a certain region. As we set the the boundary layer very small,
    # this will only be the aggregate of the granule
    aggregate, _ = determine_diffusion_region(grid2bac, grid2nBacs, bac, detachment_grid)
    biofilm = aggregate .> 0

    # Create matrices
    T = zeros(size(grid2nBacs))
    Visited = .!biofilm

    return T
end