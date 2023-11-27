function rMatrix_section(pH, conc, grid2bac, grid2nBacs, diffRegion, 
    grouped_bac, nBacs, bacOffset,
    reactive_indices, Ks, Ki, Keq, chrM, mMetabolism, mDecay, constants, kinetics)
    """
    Calculate the reaction matrix, mu and pH in a specific part of the simulation
    pH:                 Matrix of pH at each grid cell
    conc:               Matrix containing all concentrations per grid cell (ny * nx * Ncompounds)
    grid2bac:           Matrix with per gridcell which bacteria are in it
    grid2nBacs:         Matrix with per gridcell the amount of bacteria in it
    diffRegion:         Binary matrix indicating which gridcells are diffusion region
    grouped_bac:        Matrix with bacterial species, molarMass and binary activity
    nBacs:              Number of bacteria present
    bacOffset:          offset of bacteria due to chunks (=0 in sequential calculations)
    reactive_indices:   index of preferred specie
    Ks, Ki, Keq:        Ks, Ki and Keq values
    chrM:               matrix with the charge per chemical specie
    mMetabolism:        Metabolism matrix
    mDecay:             Decay matrix
    constants:          Array containing, pH_bulk, pHincluded, pHtolerance, temperature and speciation
    kinetics:           Matrix with all mu_max and maintenance

    Returns:
    reaction_matrix:    Contains per grid cell the concentrationChanges per compound
    mu:                 Vector with growth for each bacteria
    pH_new:             Matrix with updated pH
    """

    # Include files
    include(string(pwd(), "\\lib\\reaction_matrix\\solve_pH.jl"))
    include(string(pwd(), "\\lib\\reaction_matrix\\calculate_monod.jl"))
    include(string(pwd(), "\\lib\\reaction_matrix\\determine_max_growth_rate_and_maint.jl"))

    # Extract variables from parameters
    pH_bulk = constants[1]
    pHincluded = constants[2]
    pHtolerance = constants[3]
    T = constants[4]
    speciation = constants[5]

    bac_species = grouped_bac[:,1]
    bac_molarMass = grouped_bac[:,2]
    bac_active = grouped_bac[:,3]

    mu_max_list = kinetics[:,1]
    maint_list = kinetics[:,2]

    # Create storage space
    mu = zeros(nBacs)
    reaction_matrix = zeros(size(conc))
    pH_new = zeros(size(pH))

    # for each gridcell
    for x_index in axes(conc,2)
        for y_index in axes(conc, 1)
            if .!diffRegion[y_index, x_index] # If in bulk
                pH_new[y_index, x_index] = pH_bulk
                # No bacteria in bulk, so no mu or reaction_matrix update

            else # In diffusion layer, so pH calculation needs to be performed
                # Calculate pH & speciation
                if Bool(speciation)
                    Sh_old = 10^(-pH[y_index, x_index])
                    # In the following line, Concentration is a 1D vector (nComp), if problems, turn into 3D matrix (reshape(A, :, 1, 1))
                    spcM, Sh = solve_pH(Sh_old, [reshape(conc[y_index, x_index, :], :); 1; 0], Keq, chrM, pHincluded, pHtolerance) # Calculate speciation and proton concentration
                    pH_new[y_index, x_index] = -log10(Sh) # Get new pH

                else # No speciation
                    pH_new[y_index, x_index] = pH_bulk
                    Sh = 10^(-pH[y_index, x_index])
                    spcM = reshape(conc[y_index, x_index, :], :) # !1D, could be converted to a 3D matrix if necessary                
                end

                if grid2nBacs[y_index, x_index] > 0 # If cells are also found in this gridcell, update reaction matrix
                    # Get which bacteria are in this grid cell
                    iBacs = reshape(grid2bac[y_index, x_index, 1:grid2nBacs[y_index,x_index]], :) # 1D vector, could be 2D matrix

                    # Correct for chunk indexing
                    iBacs = iBacs .- bacOffset

                    speciesGrid = bac_species[iBacs] # Species for present bacteria
                    unique_species = unique(speciesGrid) # Which species are present
    
                    for curr_species in Int.(unique_species)
                        if isnan(mu_max_list[1]) # if not given as input, calculate
                            mu_max, maint = determine_max_growth_rate_and_maint(curr_species, T, Sh)
                        else
                            mu_max = mu_max_list[curr_species]
                            maint = maint_list[curr_species]
                        end
                    
                        # Get reactive concentrations for soluble components
                        reactive_conc = spcM[reactive_indices]

                        # Set mu for all bacteria of same species in that grid cell
                        M = calculate_monod(Ks[curr_species, :], Ki[curr_species, :], reactive_conc)    # Calculate the Monod-factor for mu-calculations
                        mu_noMaintenance = mu_max * M                                           # [1/h]
                        mu_withMaintenance = mu_noMaintenance - maint                           # [1/h]
                        mu[iBacs[speciesGrid .== curr_species]] .= mu_withMaintenance        # add the calculated mu at the right location to the storage

                        # Calculate cumulative mass of active bacteria in grid cell
                        cumulative_mass = sum(bac_molarMass[iBacs[speciesGrid .== curr_species]] .* bac_active[iBacs[speciesGrid .== curr_species]]) # [molX]

                        # Update reaction_matrix element for this grid cell
                        concentrationChange = mMetabolism[:, curr_species] * mu_noMaintenance           # [mol_i/molX/h] Biomass specific change
                        
                        if mu_withMaintenance < 0
                            concentrationChange = concentrationChange .- mDecay[:, curr_species] * mu_withMaintenance
                        end
                        # Total concentration Changes
                        concentrationChange = concentrationChange * cumulative_mass                     # [mol_i/h] Total change for this specie

                        # Fill in matrix
                        reaction_matrix[y_index, x_index, :] .+ reshape(concentrationChange, 1, 1, :)   # [mol_i/h]
                    end
                end
            end
        end
    end

    return reaction_matrix, mu, pH_new
end
