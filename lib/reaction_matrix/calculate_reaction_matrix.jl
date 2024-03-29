function calculate_reaction_matrix!(grid2bac, grid2nBacs, bac_vecfloat, bac_vecint, bac_vecbool, diffRegion, conc, constants_float, constants_vecfloat, constants_vecint, constants_matfloat, pH, chunks, nChunks_dir, settings_bool)
    """
    This function calculates how much of each compound is consumed per gridcell due to bacterial activity. 
    It also updates the growth rate of the respective bacteria.
    This function mainly orchestrates the creation of the reaction_matrix. It sets up the parallelization if that is desired as well.
    The actual calculation is done by rMatrix_section.jl

    Arguments
    grid2bac:           A matrix (ny, nx, ?) which contains for each gridcell which bacteria is located
                        there. The number corresponds to the index in the bac struct
    grid2nBacs:         A (ny, nx) matrix which contains for each gridcell how many bacteria are located there
    bac:                A "General" struct containing all parameters related to the bacteria
    diffRegion:         A BitMatrix indicating per gridcell whether that cell is in the diffusion region
    conc:               A (ny, nx, ncompounds) matrix containing all concentrations per gridcell
    constants:          A "General" struct containing all the simulation constants
    pH:                 A (ny, nx) matrix containing the pH value per grid cell
    chunks:             A "General" struct containing the start and end cooridinates of the chunks
    nChunks_dir:        The amount of chunks that are generated in each direction
    settings:           A "General" struct containing all the settings of the simulation

    Returns
    reaction_matrix:    A (ny, nx, ncompounds) matrix containing all reaction rates per gridcell and compound [mol/L/h]
    mu:                 An (nBacs,) vector with updated growth rates
    pH:                 A (ny, nx) matrix containing the pH value per grid cell
    """

    pHincluded = settings_bool.pHincluded
    nChunks = nChunks_dir ^2

    # convert initial pH (single value) to pH matrix
    if isa(pH, Number)
        pH = ones(size(grid2nBacs)) .* pH
    end

    # Extract variables from parameters
    Keq = constants_matfloat.Keq                                 # A (ncompounds, 4) matrix with the equilibrium constants
    chrM = constants_matfloat.chrM                               # A (ncompounds, 5) matrix with charge values
    Vg = constants_float.Vg                                   # The grid cell volume [L]
    # compoundNames = constants_vecstring.compoundNames             # A (ncompounds,) vector with the compound names (without H2O or H) COuld be needed for parallelization
    reactive_indices = constants_vecint.reactive_indices       # The indices that indicate where the reactive specie is located in the matrix
    Ks = constants_matfloat.Ks                                   # A (nSpecies, ncompounds) matrix with Ks values
    Ki = constants_matfloat.Ki                                   # A (nSpecies, ncompounds) matrix with Ki values
    mMetabolism = constants_matfloat.MatrixMet                   # A (nCompounds, nSpecies) matrix with metabolism coefficients
    mDecay = constants_matfloat.MatrixDecay                      # A (nCompounds, nSpecies) matrix with decay coefficients

    # Set up storage
    reaction_matrix = zeros(size(conc))
    mu = zeros(size(bac_vecfloat.x))

    # Pre-compute for bulk-liquid
    # At coordinates 1,1, there should never be diffusion layer
    Sh_bulk = 10^(-pH[1,1])

    if pHincluded
        # Find equilibrium between species and proton concentration for bulk and calculate pH
        _, Sh_bulk = solve_pH(Sh_bulk, [reshape(conc[1,1,:], :); 1; 0], Keq, chrM, pHincluded, constants_float.pHtolerance)
        pH_bulk = -log10(Sh_bulk)
    else
        pH_bulk = pH[1,1]
    end

    # Group constants for easy passing to several cores
    constantValues = (pH_bulk, pHincluded, constants_float.pHtolerance, constants_float.T, settings_bool.speciation)
    grouped_bac = (bac_vecint.species, bac_vecfloat.molarMass, bac_vecbool.active)
    grouped_kinetics = (constants_vecfloat.mu_max, constants_vecfloat.maintenance)

    if settings_bool.parallelized

        # ================================ PARALLEL CALCULATION START ========================================
        # Place to implement parallelization. Not implemented now as Julia might have a package that can help with this: ImplicitGlobalGrid.jl
        # ================================ PARALLEL CALCULATION END ==========================================

    else

        # ================================ SEQUENTIAL CALCULATION START ========================================
        bacOffset = 0 # Reset at 0, not used in sequential calculation
        # Calculate as a whole
        reaction_matrix, mu, pH = rMatrix_section(pH, conc, grid2bac, grid2nBacs, diffRegion, 
        grouped_bac, length(bac_vecfloat.x), bacOffset,
        reactive_indices, Ks, Ki, Keq, chrM, mMetabolism, mDecay, constantValues, grouped_kinetics)
        # ================================ SEQUENTIAL CALCULATION END ==========================================
        
    end

    # Final unit correction
    reaction_matrix = reaction_matrix ./ Vg             # [mol_i/L/h]
    return reaction_matrix, mu, pH
end
