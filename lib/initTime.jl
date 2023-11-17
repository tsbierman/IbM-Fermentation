function initTime!(grid, bac, init_params, constants, settings)
    """
    Initialises values from preset parameters for the initialisation of Time struct
    """
    include(string(pwd(), "\\lib\\Lib_Module.jl"))
    # STILL NEED TO INCLUDE THE CALCULATE_REACTION_MATRIX FILE

    # Calculate boundary conditions
    bulk_concs, invHRT = Lib_Module.calculate_bulk_concentrations(bac, constants, init_params.init_bulk_conc, init_params.invHRT, 0, constants.dT_bac, settings)

    # make bacterial-grid matrices
    grid2bac, grid2nBacs = Lib_Module.determine_where_bacteria_in_grid(grid, bac)

    # Determine diffusion layer and calculate ranges for focus mask
    diffusion_region, focus_region = Lib_Module.determine_diffusion_region(grid2bac, grid2nBacs, bac, grid)
    xRange = focus_region.x0:focus_region.x1
    yRange = focus_region.y0:focus_region.y1

    if constants.debug.plotDiffRegion
        plotDiffRegion(grid, bac, diffusion_region, true)
    end

    # Initialise concentrations and pH
    conc = zeros(grid.ny, grid.nx, length(constants.compoundNames))
    conc = Lib_Module.set_concentrations(conc, bulk_concs, diffusion_region)
    reaction_matrix = zeros(grid.ny, grid.nx, size(conc,3))
    pH = ones(grid.ny, grid.nx) .* constants.pHsetpoint

    # Set bulk layer concentrations
    conc = Lib_Module.set_concentrations(conc, bulk_concs, .!diffusion_region)

    settings.parallelized = false # at init: compute sequentially
    # Calculate reaction matrix
    # ------------Might give issues based on diffusion_region is indexed------------
    reaction_matrix[yRange,xRange, :], bac.mu, pH[yRange, xRange] = calculate_reaction_matrix!(grid2bac[yRange, xRange, :],
    grid2nBacs[yRange, xRange], bac, diffusion_region[yRange, xRange, :], conc[yRange, xRange, :], constants, pH[yRange,xRange], 0, 0, settings)

    return conc, bulk_concs, invHRT, reaction_matrix, pH, bac
end
