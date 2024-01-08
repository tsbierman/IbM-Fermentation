function initTime!(grid, bac, init_params, constants, settings)
    """
    This function initialises values from preset parameters for later initialisation of Time struct

    Arguments
    grid:               A "General" struct containing all parameters related to the grid
    bac:                A "General" struct containing all parameters related to the bacteria
    init_params:        A "General" struct containing the parameters values at the start of the simulation
    constants:          A "General" struct containing all the simulation constants
    settings:           A "General" struct containing all the settings of the simulation

    Returns
    conc:               A (ny, nx, ncompounds) matrix containing all concentrations per gridcell
    bulk_concs          A (ncompounds,) vector with the new bulk concentration
    invHRT              1 / HRT
    reaction_matrix:    A (ny, nx, ncompounds) matrix containing all reaction rates per gridcell and compound [mol/L/h]
    pH:                 A (ny, nx) matrix containing the pH value per grid cell
    bac:                A "General" struct containing all parameters related to the bacteria, updated with mu-values
    """

    # # Load required files
    # include(string(pwd(), "\\lib\\Lib_Module.jl"))
    # include(string(pwd(), "\\lib\\reaction_matrix\\calculate_reaction_matrix.jl"))

    # Calculate boundary conditions
    bulk_concs, invHRT = Lib_Module.calculate_bulk_concentrations(bac, constants, init_params.init_bulk_conc, init_params.invHRT, 0, constants.dT_bac, settings)

    # Make bacterial-grid matrices
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
    conc = Lib_Module.set_concentrations!(conc, init_params.init_concs, diffusion_region)
    reaction_matrix = zeros(grid.ny, grid.nx, size(conc,3))
    pH = ones(grid.ny, grid.nx) .* constants.pHsetpoint

    # Set bulk layer concentrations
    conc = Lib_Module.set_concentrations!(conc, bulk_concs, .!diffusion_region)

    settings.parallelized = false # at init: compute sequentially
    # Calculate reaction matrix
    reaction_matrix[yRange,xRange, :], bac.mu, pH[yRange, xRange] = calculate_reaction_matrix!(grid2bac[yRange, xRange, :],
    grid2nBacs[yRange, xRange], bac, diffusion_region[yRange, xRange, :], conc[yRange, xRange, :], constants, pH[yRange,xRange], 0, 0, settings)

    return conc, bulk_concs, invHRT, reaction_matrix, pH, bac
end
