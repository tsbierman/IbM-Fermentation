function initTime!(grid_float, grid_int, bac_vecfloat, bac_vecint, bac_vecbool, init_params, constants, settings_bool, settings_string)
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

    # Calculate boundary conditions
    bulk_concs, invHRT = calculate_bulk_concentrations(bac_vecfloat, bac_vecbool, constants, init_params.init_bulk_conc, init_params.invHRT[1], 0, constants.dT_bac, settings_bool, settings_string)

    # Make bacterial-grid matrices
    grid2bac, grid2nBacs = determine_where_bacteria_in_grid(grid_float, grid_int, bac_vecfloat)

    # Determine diffusion layer and calculate ranges for focus mask
    diffusion_region, focus_region = determine_diffusion_region(grid2bac, grid2nBacs, bac_vecfloat, grid_float, grid_int)
    xRange = focus_region.x0:focus_region.x1
    yRange = focus_region.y0:focus_region.y1

    # if constants.debug.plotDiffRegion
    #     plotDiffRegion(grid_float, bac, diffusion_region, true)
    # end

    # Initialise concentrations and pH
    conc = zeros(grid_int.ny, grid_int.nx, length(constants.compoundNames))
    conc = set_concentrations!(conc, init_params.init_concs, diffusion_region)
    reaction_matrix = zeros(grid_int.ny, grid_int.nx, size(conc,3))
    pH = ones(grid_int.ny, grid_int.nx) .* constants.pHsetpoint

    # Set bulk layer concentrations
    conc = set_concentrations!(conc, bulk_concs, .!diffusion_region)

    settings_bool.parallelized = false # at init: compute sequentially
    # Calculate reaction matrix
    reaction_matrix[yRange,xRange, :], bac_vecfloat.mu, pH[yRange, xRange] = calculate_reaction_matrix!(grid2bac[yRange, xRange, :],
    grid2nBacs[yRange, xRange], bac_vecfloat, bac_vecint, bac_vecbool, diffusion_region[yRange, xRange, :], conc[yRange, xRange, :], constants, pH[yRange,xRange], 0, 0, settings_bool)

    return conc, bulk_concs, invHRT, reaction_matrix, pH, bac_vecfloat, bac_vecint, bac_vecbool
end
