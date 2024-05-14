function integTime(simulation_file, directory)
    """
    This function integrates over Time. It regulates the whole process

    Arguments
    Simulation_file:    The file where all the variables are taken from
    Directory:          Directory where the output should be placed

    Returns
    Nothing as it saves the results in files
    """

    # Load preset file
    grid_float, grid_int, bac_vecfloat, bac_vecint, bac_vecbool, constants_float, constants_vecfloat, constants_vecint, constants_vecstring, constants_vecbool, constants_matfloat, init_params, settings_bool, settings_string = 
    load(simulation_file, "grid_float","grid_int", "bac_vecfloat", "bac_vecint", "bac_vecbool", "constants_float", "constants_vecfloat", "constants_vecint", "constants_vecstring", "constants_vecbool", "constants_matfloat", "init_params", "settings_bool", "settings_string")
    constants_vecbool.plotConvergence = [false]
    constants_vecbool.plotDiffRegion = [false]

    # -----------------PARALLELISATION --------------------------------
    # Overall settings (parallelization)
    # if settings.parallelized
    #     cores = Sys.CPU_THREADS
    #     nChunks_dir = ceil(sqrt(cores))
    #     # Checking if parpool is started or something
    #     println("Parallelization enabled for $(cores) cores\n")
    # end
    # -----------------PARALLELISATION --------------------------------

    # Determine initial simulation parameters
    backup_file = string(directory, "\\backup.jld2")
    profiling_file = string(directory, "\\profilingResults.jld2")

    if isfile(backup_file) && isfile(profiling_file)
        # load from files
        bac_vecfloat, bac_vecint, bac_vecbool, bulk_concs, invHRT, conc, reaction_matrix, pH = load(backup_file, "bac_vecfloat", "bac_vecint", "bac_vecbool", "bulk_concs", "invHRT", "conc", "reaction_matrix", "pH")
        profiling, maxErrors, normOverTime, nDiffIters, maxInitRES, bulk_history, Time, Time_vecfloat = load(profiling_file, "profiling", "maxErrors", "normOverTime", "nDiffIters", "maxInitRES", "bulk_history", "Time", "Time_vecfloat")
    else
        # Initiate from preset values
        conc, bulk_concs, invHRT, reaction_matrix, pH, bac_vecfloat, bac_vecint, bac_vecbool = initTime!(grid_float, grid_int, bac_vecfloat, bac_vecint, bac_vecbool, init_params, constants_float, constants_vecfloat, constants_vecint, constants_vecstring, constants_vecbool, constants_matfloat, settings_bool, settings_string)
        

        # Initiate time and profiling information/storage from preset
        Time = Float_struct()
        Time_vecfloat = VectorFloat_struct()
        Time.current = 0                    # Current time
        Time.steadystate = Time.current + (constants_vecint.nDiffusion_per_SScheck[1] - 1) * constants_float.dT # When to check for steadystate
        Time.save = constants_float.dT_save       # When to save
        Time.backup = constants_float.dT_backup   # When to make backup
        Time.analyse = constants_float.dT_analyse # When to analyse balances
        Time.changed_dT = 0                 # When was the latest diffusion dT change
        Time.changed_dT_bac = 0             # When was the latest bacteria dT change
        Time.dT = constants_float.dT              # Current diffusion dT
        Time.dT_bac = constants_float.dT_bac      # Current bacteria dT
        Time.bac = Time.dT_bac              # Include dT_bac and dT_divide in one variable

        if settings_bool.dynamicDT
            Time.maxDT = constants_float.maxDT
            Time.minDT = constants_float.minDT
            Time.maxDT_bac = constants_float.maxDT_bac
            Time.minDT_bac = constants_float.minDT_bac

            maximum_space_needed = ceil(Int, constants_float.simulation_end/Time.minDT_bac) + 1

            Time_vecfloat.history = zeros(Float32, maximum_space_needed)                         # Vector to save the time at each Steady-State
            profiling = zeros(Float32, maximum_space_needed, 11)                        # Matrix to save time spent on certain calculations
            maxErrors = zeros(Float32, maximum_space_needed)                            # Vector to store max Error per dT_bac
            normOverTime = zeros(Float32, maximum_space_needed)                         # Vector to store norm of concentration differance per dT_bac
            nDiffIters = zeros(UInt16, maximum_space_needed)                            # Vector to store number of diffusion iterations per steady state
            bulk_history = zeros(Float32, size(bulk_concs, 1), maximum_space_needed)    # Matrix to store bulk concentrations over time
            maxInitRES = zeros(Float32, maximum_space_needed)                           # Vector to store maximum initial RES values
        else
            max_space_needed = ceil(constants_float.simulation_end / constants_float.dT_bac) + 1

            Time_vecfloat.history = zeros(Float32, max_space_needed)                             # Vector to save the time at each Steady-State
            profiling = zeros(Float32, max_space_needed, 11)                            # Matrix to save time spent on certain calculations
            maxErrors = zeros(Float32, max_space_needed)                                # Vector to store max Error per dT_bac
            normOverTime = zeros(Float32, max_space_needed)                             # Vector to store norm of concentration differance per dT_bac
            nDiffIters = zeros(UInt16, max_space_needed)                                # Vector to store number of diffusion iterations per steady state
            bulk_history = zeros(Float32, size(bulk_concs, 1), max_space_needed)        # Matrix to store bulk concentrations over time
            maxInitRES = zeros(Float32, max_space_needed)                               # Vector to store maximum initial RES values
            Time.minDT = Time.dT
        end
        bulk_history[:,1] = bulk_concs # Is added after changing iProf, so first value should be placed already

        # Initialise saving file
        save_slice(bac_vecfloat, bac_vecint, bac_vecbool, conc, bulk_concs, pH, invHRT, 0, grid_float, grid_int, constants_float, constants_vecint, constants_vecstring, directory)
    end

    # Initialise storing space
    RESvalues = zeros(length(constants_vecstring.compoundNames), 10000) # Reserve space for n steady state checks beforehand (can be more)
    norm_diff = zeros(10000)
    res_bacsim = zeros(10000, 2)

    iProf = findfirst(profiling .== 0)[1]       # Keep track of index of profiling (every simulated dT_bac +1 index) (starts half way if restarting from storage)
    iDiffusion = 1                              # Keep track of index of diffusion (per 1 dT_bac: iDiffusion == cycles of diffusion)
    iRES = 0                                    # Times steady state has been calculated this dT_bac

    # Make bacterial grid-matrix
    grid2bac, grid2nBacs = determine_where_bacteria_in_grid(grid_float, grid_int, bac_vecfloat)

    # Determine diffusion layer and calculate ranges for focus mask
    diffusion_region, focus_region = determine_diffusion_region(grid2bac, grid2nBacs, bac_vecfloat, grid_float, grid_int)
    xRange = focus_region.x0:focus_region.x1
    yRange = focus_region.y0:focus_region.y1

    # ----------------- START PARALLELISATION --------------------------------
    if settings_bool.parallelized

        # Create chunks
        chunks_matrix, chunks_int = create_chunks(nChunks_dir, focus_region)

        # sort bacteria
        bac_vecfloat, bac_vecint, bac_vecbool = sort_bacteria_into_chunks!(bac_vecfloat, bac_vecint, bac_vecbool, grid_float, chunks_int, focus_region, nChunks_dir)

        # Recalculate the grid2bac matrix
        grid2bac, _ = determine_where_bacteria_in_grid(grid_float, grid_int, bac_vecfloat)
    else
        # set dummy values for chunk variables
        chunks = 0
        nChunks_dir = 0
    end
    # ------------------- END PARALLELISATION --------------------------------

    # Time advancements (dT / dT_steadystate)
    prev_conc = copy(conc) # Store current concs before calculating diffusion

    while Time.current < constants_float.simulation_end

        if mod(iDiffusion, constants_vecint.nItersCycle[1]) == 0
            println("Currently at diffusion iteration $(iDiffusion) (max error: $(maximum(RESvalues[:, iRES])))\n")
        end

        # diffuse (MG)
        tick()
        try
            conc[yRange, xRange, :] = diffusionMG!(conc[yRange, xRange, :], reaction_matrix[yRange, xRange, :], bulk_concs, diffusion_region[yRange, xRange], grid_float, constants_float, constants_vecfloat, constants_vecstring, Time)
        catch e
            if isa(e, ErrorException)
                Time = decrease_dT_diffusion!(Time, getfield(e, :msg), grid_float.dx, constants_vecfloat)
            else
                rethrow(e)
            end
        end
        profiling[iProf, 1] = profiling[iProf, 1] + tok()

        # Set bulk layer concentrations (in theory, not needed anymore with correct diffusion model)
        conc = set_concentrations!(conc, bulk_concs, .!diffusion_region)

        # Update bacterial mass
        tick()
        bac_vecfloat = update_bacterial_mass!(bac_vecfloat, bac_vecbool, Time.dT)
        profiling[iProf, 2] = profiling[iProf, 2] + tok()

        # Calculate reaction_matrix
        tick()
        reaction_matrix[yRange, xRange, :], bac_vecfloat.mu, pH[yRange, xRange] = calculate_reaction_matrix!(grid2bac[yRange, xRange, :],
        grid2nBacs[yRange, xRange], bac_vecfloat, bac_vecint, bac_vecbool, diffusion_region[yRange, xRange, :], conc[yRange, xRange, :], 
        constants_float, constants_vecfloat, constants_vecint, constants_matfloat, pH[yRange, xRange], chunks, nChunks_dir, settings_bool)
        profiling[iProf, 3] = profiling[iProf, 3] + tok()

        # If T > T_ss: calcualte residual
        if Time.current >= Time.steadystate

            # Increase counter of RES checks performed
            iRES = iRES + 1

            # Perform check for steady state
            tick()
            ssReached, RESvalues[:, iRES] = steadystate_is_reached(conc[yRange, xRange, :], reaction_matrix[yRange, xRange, :], grid_float.dx, bulk_concs, diffusion_region[yRange, xRange], constants_float, constants_vecfloat, constants_vecstring)
            norm_diff[iRES] = sqrt(sum((prev_conc .- conc).^2)) # Difference between two diffusions
            res_bacsim[iRES, 1] = maximum(abs.((prev_conc .- conc) ./ Time.dT)) # maximum difference relative to passed time
            res_bacsim[iRES, 2] = norm_diff[iRES] / Time.dT                     # Normalised difference relative to passed time
            profiling[iProf, 4] = profiling[iProf, 4] + tok()
            
            prev_conc = copy(conc)  # Store the newly calculated conc before changing again in next cycle 
            
            # Perform dynamic dT for diffusion
            if settings_bool.dynamicDT
                if slow_convergence(iRES, RESvalues, constants_float, constants_vecint) && Time.dT < Time.maxDT/2
                    # Time = increase_dT_diffusion!(Time, "Slow convergence", grid.dx, constants_vecfloat)
                    Time = decrease_dT_diffusion!(Time, "Diffusion takes a long time", grid_float.dx, constants_vecfloat)
                end

                if upward_trend(iRES, RESvalues)
                    Time = decrease_dT_diffusion!(Time, "Upward trend in RES values detected", grid_float.dx, constants_vecfloat)
                    Time = decrease_dT_diffusion!(Time, "Upward trend in RES values detected", grid_float.dx, constants_vecfloat)
                elseif non_convergent(iRES, RESvalues, constants_float.tolerance_no_convergence)
                    Time = decrease_dT_diffusion!(Time, "Convergence is stuck", grid_float.dx, constants_vecfloat)
                end
            else
                if iDiffusion > 5000 && slow_convergence(iRES, RESvalues, constants_float, constants_vecint)
                    # Without dynamic timestep & negative concentrations
                    # Due to too large step size, accept SS under non-convergent conditions
                    ssReached = true
                end

                if iDiffusion > 10000
                    ssReached = true
                end
            end

            if ssReached
                @printf("Steady state reached after %d diffusion iteration\n", iDiffusion)
                @printf("\twith at most %.4g mol/L/h off of steady state (norm = %e)\n", maximum(RESvalues[:,iRES]), norm_diff[iRES])

                # Set time to next bacterial activity time (skipping ahead as we have steady state)
                previousTime = Time.current         # Save time where we found steadystate
                Time.current = Time.bac

                if Time.current > constants_float.simulation_end
                    Time.current = constants_float.simulation_end - Time.dT # Prevent another bacterial timestep
                end

                # Perform dynamic dT for diffusion (for next iteration)
                if settings_bool.dynamicDT && multiple_high_iters(iDiffusion, iProf, nDiffIters, Time, constants_vecint)
                    Time = increase_dT_diffusion!(Time, "Multiple steady states reached with more than $(constants_vecint.iterThresholdIncrease[1]) diffusion iterations", grid_float.dx, constants_vecfloat)
                end

                # Calculate actual dT for integration of bacterial mass (for skipped time)
                dT_actual = Time.current - previousTime

                # Time advancements (dT_bac)
                if Time.current >= Time.bac
                    # Store values at time of reaching steady state
                    maxErrors[iProf] = maximum(RESvalues[:, iRES])
                    normOverTime[iProf] = norm_diff[iRES]
                    nDiffIters[iProf] = iDiffusion
                    maxInitRES[iProf] = maximum(RESvalues[:,1])
                    # Reset counters for next iteration
                    iDiffusion = 1
                    iRES = 0
                    RESvalues = zeros(length(constants_vecstring.compoundNames), 10000)

                    # Reaction_matrix & mu & pH are already calculated (steady state so still valid)

                    # Update bacteria: mass
                    tick()
                    bac_vecfloat = update_bacterial_mass!(bac_vecfloat, bac_vecbool, dT_actual)
                    profiling[iProf, 2] = profiling[iProf, 2] + tok()

                    # Time advancements (dT_divide)
                    # Determine radius bacteria from mass
                    tick()
                    bac_vecfloat = update_bacterial_radius!(bac_vecfloat, constants_float)
                    profiling[iProf, 2] = profiling[iProf, 2] + tok()

                    # Bacteria: inactivate or die
                    tick()
                    if constants_vecbool.inactivationEnabled[1]
                        bac_vecbool = bacteria_inactivate!(bac_vecfloat, bac_vecbool, constants_float)
                    else
                        bac_vecfloat, bac_vecint, bac_vecbool = bacteria_die!(bac_vecfloat, bac_vecint, bac_vecbool, constants_float)
                    end

                    # Bacteria: Divide
                    bac_vecfloat, bac_vecint, bac_vecbool, nDivCycles = bacteria_divide!(bac_vecfloat, bac_vecint, bac_vecbool, constants_float)
                    if nDivCycles > 1
                        Time = decrease_dT_bac!(Time, "Bacteria are dividing too fast")
                    end
                    profiling[iProf, 5] = profiling[iProf, 5] + tok()

                    # Shove bacteria
                    tick()
                    bac_vecfloat = bacteria_shove!(bac_vecfloat, grid_float, grid_int, constants_float)
                    bac_vecfloat = bacteria_shove!(bac_vecfloat, grid_float, grid_int, constants_float) # Second shove in case of overcrowding
                    profiling[iProf, 6] = profiling[iProf, 6] + tok()
                    
                    # Bacteria: detachment (for now only rough detachment is implemented)
                    tick()
                    bac_vecfloat, bac_vecint, bac_vecbool = bacteria_detachment!(bac_vecfloat, bac_vecint, bac_vecbool, grid_float, grid_int, constants_float, settings_string, Time.dT_bac, invHRT)
                    profiling[iProf, 7] = profiling[iProf, 7] + tok()

                    # Display number of bacteria in system
                    println("Current number of bacteria $(length(bac_vecfloat.x)) ($(sum(bac_vecbool.active)) active)")

                    # Auto detect when to switch to parallel compution of the rMatrix
                    # TODO: determine exact cutoff value (will be around 15000 approx)
                    if settings_bool.parallelized == false && length(bac_vecfloat.x) > 15000
                        cores = Sys.CPU_THREADS
                        nChunks_dir = ceil(sqrt(cores))
                        # DO PARALLELISATION initiation
                        println("Parallelization enabled for $(cores) cores\n")
                    end

                    # Update/re-determine where bacs are (as we have had division and shoving)
                    tick()
                    grid2bac, grid2nBacs = determine_where_bacteria_in_grid(grid_float, grid_int, bac_vecfloat)
                    profiling[iProf, 8] = profiling[iProf, 8] + tok()

                    # Update diffusion region (also changed)
                    tick()
                    diffusion_region, focus_region = determine_diffusion_region(grid2bac, grid2nBacs, bac_vecfloat, grid_float, grid_int)
                    xRange = focus_region.x0:focus_region.x1
                    yRange = focus_region.y0:focus_region.y1
                    profiling[iProf, 9] = profiling[iProf, 9] + tok()

                    if settings_bool.parallelized
                        tick()
                        # create chunks
                        chunks_matrix, chunks_int = create_chunks(nChunks_dir, focus_region)

                        # sort bacteria
                        bac_vecfloat, bac_vecint, bac_vecbool = sort_bacteria_into_chunks!(bac_vecfloat, bac_vecint, bac_vecbool, grid_float, chunks_int, focus_region, nChunks_dir)
                        profiling[iProf, 11] = profiling[iProf, 11] + tok()

                        # recalculate the grid2bac matrix
                        grid2bac, _ = determine_where_bacteria_in_grid(grid_float, grid_int, bac_vecfloat)
                    end

                    # Apply dynamic dT_bac
                    if settings_bool.dynamicDT && multiple_low_initRES(iProf, maxInitRES, Time, constants_float, constants_vecint)
                        Time = increase_dT_bac!(Time, "Multiple steady state cycles with initRES value below $(round(constants_float.initRESThresholdIncrease * 100, digits=1))")
                    end

                    # Prepare for next steadystate cycle
                    println("\n==============================\n")
                    println("Current simulation time: $(Time.current) h \n\n")

                    # Recompute reaction matrix for the next cycle
                    tick()
                    reaction_matrix[yRange, xRange, :], bac_vecfloat.mu, pH[yRange, xRange] = calculate_reaction_matrix!(grid2bac[yRange, xRange, :],
                    grid2nBacs[yRange, xRange], bac_vecfloat, bac_vecint, bac_vecbool, diffusion_region[yRange, xRange, :], conc[yRange, xRange, :], 
                    constants_float, constants_vecfloat, constants_vecint, constants_matfloat, pH[yRange, xRange], chunks, nChunks_dir, settings_bool)
                    profiling[iProf, 3] = profiling[iProf, 3] + tok()

                    # Calculate and set bulk concentrations
                    tick()
                    new_bulk_concs = copy(prev_conc)
                    while true
                        new_bulk_concs, invHRT = calculate_bulk_concentrations(bac_vecfloat, bac_vecbool, constants_float, constants_vecfloat, constants_vecint, constants_vecstring, constants_vecbool, constants_matfloat, bulk_concs, invHRT, reaction_matrix, Time.dT_bac, settings_bool, settings_string)
                        if !settings_bool.dynamicDT || bulk_conc_diff_within_limit(new_bulk_concs, bulk_concs, constants_float)
                            break
                        end

                        if Time.dT_bac <= Time.minDT_bac
                            println("Smallest dT_bac reached, cannot decrease smaller than $(Time.minDT_bac) h\n")
                            break
                        end

                        Time = decrease_dT_bac!(Time, "Too large bulk concentration jump detected")
                    end

                    bulk_change = (new_bulk_concs .- bulk_concs) ./ Time.dT_bac # [mol_i/L/h]

                    bulk_concs = copy(new_bulk_concs)
                    conc = set_concentrations!(conc, bulk_concs, .!diffusion_region)
                    profiling[iProf, 10] = profiling[iProf, 10] + tok()

                    # Place for balance check/analyse metabolites
                    if Time.current >= Time.analyse

                        biomass_close, balance_close, dirichlet_close = check_balances(bac_vecfloat, bac_vecint, bac_vecbool, constants_float, constants_vecfloat, constants_vecint, constants_vecstring, constants_vecbool, constants_matfloat, settings_string, reaction_matrix, bulk_concs, invHRT, bulk_change, 1e-3)
                        println("Biomass closes: $(biomass_close)")
                        println("Balances closing: $(balance_close)")
                        println("Dirichlet closing: $(dirichlet_close)")

                        Time.analyse = Time.analyse + constants_float.dT_analyse
                    end

                    iProf = iProf + 1
                    bulk_history[:, iProf] = bulk_concs
                    Time_vecfloat.history[iProf] = Time.current

                    # Set next bacterial time
                    Time.bac = Time.bac + Time.dT_bac
                end

                # Time advancements (dT_save)
                if Time.current >= Time.save
                    # Set next save time
                    Time.save = Time.save + constants_float.dT_save

                    # Save all important variables
                    # save_slice(bac_vecfloat, bac_vecint, bac_vecbool, conc, bulk_concs, pH, invHRT, Time.current, grid_float, grid_int, constants_float, constants_vecint, constants_vecstring, directory)
                    save_profile(bac_vecfloat, bac_vecint, bac_vecbool, conc, bulk_concs, pH, invHRT, Time.current, grid_float, grid_int, constants_float, constants_vecint, constants_vecstring, directory) # Entire plane of simulation

                    if Time.current >= Time.backup
                        # Set next backup time
                        Time.backup = Time.backup + constants_float.dT_backup

                        # Save all important variables for continuing simulation from this profilingResults
                        save_backup(bac_vecfloat, bac_vecint, bac_vecbool, bulk_concs, invHRT, conc, reaction_matrix, pH, directory)
                        save_profiling(profiling, maxErrors, normOverTime, nDiffIters, maxInitRES, bulk_history, Time, Time_vecfloat, directory)
                    end

                    if settings_string.detachment == "SBR"
                        bac_vecfloat.molarMass[bac_vecbool.active] = maximum([bac_vecfloat.molarMass[bac_vecbool.active] .- 0.1*constants_float.min_bac_mass_grams/constants_float.bac_MW, eps])
                    end
                end
            end
            # Set next steadystate time
            Time.steadystate = Time.current + constants_vecint.nDiffusion_per_SScheck[1] * Time.dT
        end

        # Post-dT updates
        # Advance current simulation time
        Time.current = Time.current + Time.dT
        iDiffusion = iDiffusion + 1
    end

    # Save all important variables one last time?
    # save_slice(bac_vecfloat, bac_vecint, bac_vecbool, conc, bulk_concs, pH, invHRT, Time.current, grid_float, grid_int, constants_float, constants_vecint, constants_vecstring, directory)         # Slice of simulation
    save_profile(bac_vecfloat, bac_vecint, bac_vecbool, conc, bulk_concs, pH, invHRT, Time.current, grid_float, grid_int, constants_float, constants_vecint, constants_vecstring, directory)       # Entire plane of simulation
    save_backup(bac_vecfloat, bac_vecint, bac_vecbool, bulk_concs, invHRT, conc, reaction_matrix, pH, directory)                      # Backup to start up halfway
    save_profiling(profiling, maxErrors, normOverTime, nDiffIters, maxInitRES, bulk_history, Time, Time_vecfloat, directory)   # Save performance
end
