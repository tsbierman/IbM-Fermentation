function determine_where_bacteria_in_grid(grid, bac)
    """
    This functions determines which gridcel contain bacteria and which bacteria
    are located in each gridcell

    Arguments
    grid:               A "General" struct containing all parameters related to the grid
    bac:                A "General" struct containing all parameters related to the bacteria

    Returns
    grid2bac:           A matrix (ny, nx, ?) which contains for each gridcell which bacteria is located
                        there. The number corresponds to the index in the bac struct
    grid2nBacs:         A (ny, nx) matrix which contains for each gridcell how many bacteria are located there
    """
    maxBacPerGrid = 4 # Based on radius, this number is maximum 19 (all minimum radius)

    # Determine grid cell per bacteria
    index_x = ceil.(bac.x / grid.dx)
    index_y = ceil.(bac.y / grid.dy)
    bac_grid = hcat(index_x, index_y)

    # Generate matrices
    grid2nBacs = zeros(UInt16, grid.ny, grid.nx) # UInt8 for reduced storage requirements (previously exceeded 255 bac per cell, therefor UInt8 not enough)
    grid2bac = zeros(UInt32, grid.ny, grid.nx, maxBacPerGrid) # UInt32 for reduced storage requirements --> lower storage requirements than sparse matrix (nBacs * nGridcells)

    # Fill matrices
    for i in axes(bac_grid)[1]
        ix = Int(bac_grid[i, 1])
        iy = Int(bac_grid[i, 2])
        grid2bac[iy, ix, grid2nBacs[iy, ix] + 1] = i
        grid2nBacs[iy, ix] = grid2nBacs[iy, ix] + 1
    end

    return grid2bac, grid2nBacs
end
