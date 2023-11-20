function determine_where_bacteria_in_grid(grid, bac)
    """
    This functions generates two matrices:
    grid2bac is a nx * ny * ? matrix with a reference which bacteria reside in which grid cell
    Thus, the values are indices of the bacteria in the bac struct
    grid2nBacs is a nx * ny matrix which contains the number of bacteria in each of the grid cells

    grid is a struct containing all static information about the grid
    bac is a struct containing all information regarding the bacteria
    """
    maxBacPerGrid = 4 # Based on radius, this number is maximum 19 (all minimum radius)

    # Determine grid cell per bacteria
    index_x = ceil.(bac.x / grid.dx)
    index_y = ceil.(bac.y / grid.dy)
    bac_grid = hcat(index_x, index_y)

    # Generate matrices
    grid2nBacs = zeros(UInt16, grid.ny, grid.nx) # UInt8 for reduced storage requirements (previously exceeded 255 bac per cell, therefor uint8 not enough)
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
