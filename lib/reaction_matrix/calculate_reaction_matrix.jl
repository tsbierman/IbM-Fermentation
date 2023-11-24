function calculate_reaction_matrix(grid2bac, grid2nBacs, bac, diffRegion, conc, constants, pH, chunks, nChunks_dir, settings)
    """
    Calculate how many compoundNames are consumed per grid cell due to bacterial activity. 
    Also updates the growth rate of the respective bacteria.
    This function mainly orchestrates the creation of the reaction_matrix. It sets up the parallelization is desired as well.
    The actual calculation is done by rMatrix_section.jl

    grid2bac:           matrix with per grid cell which bacteria reside there
    grid2nBacs:         matrix with how many bacteria per grid cell
    bac:                struct containing all information regarding the bacteria
    diffRegion:         The coordinates of the region that has to be calculated
    conc:               matrix containing all concentrations per grid cell as of (ix, iy, compound)
    constants:          struct containing all simulation constants
    pH:                 matrix with per grid cell the pH (or single value for whole region)
    chunks:             struct containing info regarding the chunks (for parallelization)
    nChunks_dir:        the number of chunks per direction (for parallelization)
    settings:           struct containing simulation wide settings

    Returns:
    reaction_matrix:    matrix with per grid cell and per compound the change [h-1] due to bacterial activity
    mu:                 vector with updates growth rates per bacteria
    pH:                 matrix with per grid cell the pH
    """

    # Including files
    include(string(pwd(), "\\lib\\reaction_matrix\\solve_pH.jl"))
    include(string(pwd(), "\\lib\\reaction_matrix\\rMatrix_section.jl"))

    pHincluded = settings.pHincluded
    nChunks == nChunks_dir ^2

    # convert init pH (single value) to pH matrix
    if isscalar(pH)
        pH = ones(size(grid2nBacs)) * pH
    end

    # Extract variables from parameters
    Keq = constants.Keq                                 # Equilibrium constant matrix
    chrM = constants.chrM                               # Charge matrix
    Vg = constants.Vg                                   # Grid cell volume [L]
    compoundNames = constants.compoundNames             # Compound names (without H2O or H)
    reactive_indices = constants.reactive_indices       # Indices that indicate where the reactive specie is located in the matrix
    Ks = constants.Ks                                   # Matrix with Ks values [nBacSpecies * nCompounds]
    Ki = constants.Ki                                   # Matrix with Ki values [nBacSpecies * nCompounds]
    mMetabolism = constants.MatrixMet                   # Matrix with metabolism coefficients [nCompounds * nBacSpecies]
    mDecay = constants.MatrixDecay                      # Matrix with decay coefficients [nCompounds * nBacSpecies]

    # Set up storage
    reaction_matrix = zeros(size(conc))
    mu = zeros(size(bac.x))

    # Pre-compute for bulk-liquid
    # At coordinates 1,1, there should never be diffusion layer
    Sh_bulk = 10^(-pH[1,1])

    if pHincluded
        # Find equilibrium between species and proton concentration for bulk and calculate pH
        _, Sh_bulk = solve_pH(Sh_bulk, [reshape(conc[1,1,:], :); 1; 0], Keq, chrM, pHincluded, constants.pHtolerance)
        pH_bulk = -log10(Sh_bulk)
    else
        pH_bulk = pH[1,1]
    end

    # Group constants for easy passing to several cores
    constantValues = [pH_bulk, pHincluded, constants.pHtolerance, constants.T, settings.speciation]
    grouped_bac = [bac.species bac.molarMass bac.active]
    grouped_kinetics = [constants.mu_max constants.maintenance]

    if settings.parallelized
        # ================================ PARALLEL CALCULATION START ========================================
        # Place to implement parallelization. Not implemented now as Julia might have a package that can do this: ImplicitGlobalGrid.jl
        # ================================ PARALLEL CALCULATION END ========================================
    else
        # ================================ SEQUENTIAL CALCULATION START ========================================
        bacOffset = 0 # Reset at 0, not used in sequential calculation
        # Calculate as a whole
        reaction_matrix, mu, pH = rMatrix_section(pH, conc, grid2bac, grid2nBacs, diffRegion, 
        grouped_bac, length(bac.x), bacOffset,
        reactive_indices, Ks, Ki, Keq, chrM, mMetabolism, mDecay, constantValues, grouped_kinetics)
        # ================================ SEQUENTIAL CALCULATION END ========================================
    end

    # Final unit correction
    reaction_matrix = reaction_matrix ./ Vg
    return reaction_matrix, mu, pH
end
