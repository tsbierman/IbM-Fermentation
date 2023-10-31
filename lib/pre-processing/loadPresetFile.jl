function loadPresetFile(filename)
    """
    Reads the excel and stores the information in structs for later usage
    """
    # Initialisation of the structs
    grid = General()
    constants = General()
    settings = General()
    init_params = General()
    bac_init = General()

    #Read file
    file = XLSX.readxlsx(filename)

    # initialise grid
    names_discr, values_discr = collect(skipmissing(file["Discretization"][:,1])), collect(skipmissing(file["Discretization"][:,2]))
    grid.dx = values_discr[names_discr .== "dx"][1]                                        # [m]
    grid.dy = values_discr[names_discr .== "dy"][1]                                        # [m]
    grid.nx = values_discr[names_discr .== "nx"][1]                                        # [-]
    grid.ny = values_discr[names_discr .== "ny"][1]                                        # [-]
    grid.blayer_thickness = values_discr[names_discr .== "Boundary layer thickness"][1]    # [m]

    constants.Vg = (grid.dx ^ 3) * 1000                                     # [L]
    constants.max_granule_radius = ((grid.nx - 4) * grid.dx) / 2            # [m]

    # initialise constants (Time)
    constants.simulation_end = values_discr[names_discr .== "Simulation end"][1]           # [h]
    constants.dT = values_discr[names_discr .== "Initial dT diffusion"][1]                 # [h]
    constants.dT_bac = values_discr[names_discr .== "Initial dT bacteria"][1]              # [h]
    constants.dT_save = values_discr[names_discr .== "dT save"][1]                         # [h]
    constants.dT_backup = values_discr[names_discr .== "dT backup"][1]                     # [h]

    settings.dynamicDT = values_discr[names_discr .== "Dynamic dT"][1]                  # [Bool]

    # Only if dynamic time stepping is enabled
    if settings.dynamicDT
        dynamicDT = General()   # Extra structure to store in the other structure
        dynamicDT.nIterThreshold = values_discr[names_discr .== "nIterThreshold"][1]
        dynamicDT.iterThresholdDecrease = values_discr[names_discr .== "iterThresholdDecrease"][1]
        dynamicDT.iterThresholdIncrease = values_discr[names_discr .== "iterThresholdIncrease"][1]
        dynamicDT.initRESThresholdIncrease = values_discr[names_discr .== "initial RES threshold increase"][1]
        dynamicDT.nItersCycle = values_discr[names_discr .== "nIters per cycle"][1]
        dynamicDT.tolerance_no_convergence = values_discr[names_discr .== "tolerance no-convergence"][1]
        dynamicDT.maxRelDiffBulkConc = values_discr[names_discr .== "maximum relative bulk conc change"][1]
        
        dynamicDT.maxDT = values_discr[names_discr .== "Maximum dT diffusion"][1]          # [h]
        dynamicDT.minDT = values_discr[names_discr .== "Minimum dT diffusion"][1]          # [h]
        dynamicDT.maxDT_bac = values_discr[names_discr .== "Maximum dT bacteria"][1]       # [h]
        dynamicDT.minDT_bac = values_discr[names_discr .== "Minimum dT bacteria"][1]       # [h]

        constants.dynamicDT = dynamicDT
    end

    # Constants (Diffusion)
    names_diff, values_diff = collect(skipmissing(file["Diffusion"][:,1])), collect(skipmissing(file["Diffusion"][:,2]))
    constants.compoundNames = names_diff
    nCompounds = length(constants.compoundNames)
    constants.diffusion_rates = values_diff                                      # [m2/h]

    # Constants (Operational parameters)
    names_para, values_para = collect(skipmissing(file["Parameters"][:,1])), collect(skipmissing(file["Parameters"][:,2])) # It takes some extra empty rows, this removes that
    constants.pHsetpoint = values_para[names_para .== "pH setpoint"][1]                  # [-]
    constants.T = values_para[names_para .== "Temperature (K)"][1]                       # [K]
    constants.Vr = values_para[names_para .== "Representative volume"][1] * 1000         # [L]
    constants.reactor_density = values_para[names_para .== "density reactor"][1]         # [g/L]

    settings.variableHRT = values_para[names_para .== "Variable HRT"][1]                 # [Bool]
    init_params.invHRT = 1 / values_para[names_para .== "HRT"][1]                        # [1/h]

    # Only if varaible HRT is enabled
    if settings.variableHRT
        constants.bulk_setpoint = values_para[names_para .== "Setpoint"][1]              # [mol/L]
        compound_name = values_para[names_para .== "Compound setpoint"][1] 
        constants.setpoint_index = findall(constants.compoundNames .== compound_name) # [CartasianIndex]
    end

    # Constants (Bacteria)
    names_bac, values_bac = collect(skipmissing(file["Bacteria"][:,1])), collect(skipmissing(file["Bacteria"][:,2]))  # It takes some extra empty rows, this removes that
    constants.bac_MW = values_bac[names_bac .== "Molecular weight bacterium"][1]       # [g/mol]
    constants.bac_rho = values_bac[names_bac .== "Density bacterium"][1]               # [g/m3]
    constants.max_nBac = values_bac[names_bac .== "Maximum nBacteria"][1]              # [-] Maybe find better way of calculating this
    constants.inactivationEnabled = values_bac[names_bac .== "Inactivation enabled"][1] # [Bool]
    constants.min_bac_mass_grams = values_bac[names_bac .== "Minimum mass bacterium"][1] # [g]
    constants.max_bac_mass_grams = values_bac[names_bac .== "Maximum mass bacterium"][1] # [g]
    constants.bac_max_radius = values_bac[names_bac .== "Maximum radius bacterium"][1] # [m]
    constants.kDist = values_bac[names_bac .== "kDist"][1]                             # [-]
    constants.max_granule_radius = min(constants.max_granule_radius, values_bac[names_bac .== "Maximum granule radius"][1]) # [m] Either based on grid size or set value
    constants.kDet = values_bac[names_bac .== "Detachment constant"][1]                # [1/m2.h]
    settings.detachment = values_bac[names_bac .== "Detachment method"][1]          # [mechanistic, naive, none]


    # Constants (Solver)
    names_solv, values_solv = collect(skipmissing(file["Solver"][:,1])), collect(skipmissing(file["Solver"][:,2]))
    constants.diffusion_accuracy = values_solv[names_solv .== "Diffusion tolerance"][1]              # [-]
    constants.steadystate_tolerance = values_solv[names_solv .== "Steady state RES threshold"][1]    # [mol/L/h]
    tol_abs = values_solv[names_solv .== "Concentration tolerance"][1]                               # [mol/L]
    constants.correction_concentration_steady_state = tol_abs / constants.steadystate_tolerance
    constants.RESmethod = values_solv[names_solv .== "RES determination method"][1]

    constants.nDiffusion_per_SScheck = values_solv[names_solv .== "nIters diffusion per SS check"][1]

    settings.pHbulkCorrection = values_solv[names_solv .== "pH bulk concentration corrected"][1]
    settings.pHincluded = values_solv[names_solv .== "pH solving included"][1]
    # Speciation is always included if pH is included, else, read Boolean from excel
    settings.speciation = values_solv[names_solv .== "Speciation included"][1] || settings.pHincluded

    if settings.pHincluded
        constants.pHtolerance = values_solv[names_solv .== "pH solver tolerance"][1]
    else
        constants.pHtolerance = NaN
    end

    settings.structure_model = values_solv[names_solv .== "Structure model"][1]
    if settings.structure_model
        settings.type = values_solv[names_solv .== "Structure model type"][1]
    end

    settings.parallelized = values_solv[names_solv .== "Parallelisation"][1]

    # Constants (influent)
    names_temp, values_temp, condition_type_temp = file["Influent"][:,1], file["Influent"][:,2], file["Influent"][:,4] 
    names_infl, values_infl, condition_type = collect(skipmissing(names_temp)), collect(skipmissing(values_temp)), collect(skipmissing(condition_type_temp)) # It takes some extra empty rows, this removes that

    if names_infl != constants.compoundNames
        throw(ErrorException("Compounds do not have the same name or are not in the same order"))
    end

    constants.Dir_k = condition_type .== "D"
    constants.influent_concentrations = values_infl                                  # [mol/L]

    # Constants (initial conditions)
    names_init, values_init = collect(skipmissing(file["Initial condition"][:,1])), collect(skipmissing(file["Initial condition"][:,2]))

    if names_init != constants.compoundNames
        throw(ErrorException("Compounds do not have the same name or are not in the same order"))
    end

    # For the species with a Dirichlet, the influent concentration is taken as the bulk concentration
    init_params.init_concs = values_init                                             # [mol/L]
    init_params.init_bulk_conc = init_params.init_concs                         # [mol/L]
    init_params.init_bulk_conc[constants.Dir_k] = constants.influent_concentrations[constants.Dir_k]

    # Constants (Equilibrium constants & charge matrix)
    thermodynamic_parameters = file["ThermoParam"]
    names_thermo = thermodynamic_parameters[:,1]
    nColumns = length(file["ThermoParam"][2,:])-4 # -4 due to Compound names, preferred subspecies, compound phase and cell with text 

    # Divide in sections and locate H2O and H indices
    section_starts = findall(names_thermo .== constants.compoundNames[1])
    section_ends = findall(names_thermo .== constants.compoundNames[end])
    H2O_index = findall(names_thermo .== "H2O")
    H_index = findall(names_thermo .== "H")
    dG_rows = collect(section_starts[1]:section_ends[1])

    # Check the order of the compounds names
    compound_names_matrix = reshape(constants.compoundNames, length(constants.compoundNames),:) # Change shape to compare with thermodynamic names
    if thermodynamic_parameters[:,1][dG_rows] != (compound_names_matrix)
        throw(ErrorException("Compounds do not have the same name or are not in the same order"))
    end

    # Combine a section with the H2O and H indices
    dG_rows = push!(vec(dG_rows), H2O_index[1], H_index[1])

    preferred_state = thermodynamic_parameters[:,8][dG_rows]

    # Get next section indices
    Keq_start = section_starts[2][1]
    Keq_end = section_ends[2][1] + 2 # For adding H20 and H    
    constants.Keq = thermodynamic_parameters[Keq_start:Keq_end, 2:5]

    # And the last section indices
    charge_start = section_starts[3][1]
    charge_end = section_ends[3][1] 
    charge_rows = collect(section_starts[3]:section_ends[3])

    # Check order of compounds
    if thermodynamic_parameters[:,1][charge_rows] != compound_names_matrix
        throw(ErrorException("Compounds do not have the same name or are not in the same order"))
    end

    # Extract charge matrix and substitute everything without charge
    charge_end = charge_end + 2 # For adding H2O and H
    constants.chrM = thermodynamic_parameters[charge_start:charge_end, 2:6]
    idx = findall(constants.chrM .== "NA")
    constants.chrM[idx] .= 0

    # Constants (Ks and Ki)
    ks_file = file["Ks"]
    ki_file = file["Ki"]

    # Extract and check microbial species Names
    constants.speciesNames = collect(skipmissing(ks_file[:,1][2:end,1]))

    if collect(skipmissing(ki_file[:,1][2:end,1])) != constants.speciesNames
        throw(ErrorException("Bacteria species do not have the same name or are not in the same order"))
    end

    # Get names of compounds, make a unique set and extract matrix of values
    comp_names_ks = collect(skipmissing(ks_file[1,1:end-1]))
    comp_names_ki = collect(skipmissing(ki_file[1,1:end-1]))
    uniq_compounds = unique([comp_names_ks; comp_names_ki])
    values_ks = reshape(collect(skipmissing(ks_file[2:end,2:end-1])),length(constants.speciesNames), length(comp_names_ks))
    values_ki = reshape(collect(skipmissing(ki_file[2:end,2:end-1])),length(constants.speciesNames), length(comp_names_ki))

    constants.Ks = zeros(length(constants.speciesNames),length(uniq_compounds))
    constants.Ki = zeros(length(constants.speciesNames),length(uniq_compounds))

    # Loop to check whether a compound has a Ks and Ki, if yes, store value
    for (index, uniq) in enumerate(uniq_compounds)
        columnidx = findall(uniq .== comp_names_ks)
        if length(columnidx) != 0
            constants.Ks[:,columnidx] = values_ks[:,columnidx]
        end

        columnidx = findall(uniq .== comp_names_ki)
        if length(columnidx) != 0
            constants.Ki[:,columnidx] = values_ki[:,columnidx]
        end
    end

    # Reactive indices
    constants.reactive_indices = zeros(length(uniq_compounds), 1)
    # TODO Need to consider whether the LinearIndices are actually required or that it could also be done with CartesianIndex
    if settings.speciation
        sz = (nCompounds+2, nColumns)

        for (index, uniq) in enumerate(uniq_compounds)
            I = findall(constants.compoundNames .== uniq)[1][1]
            J = preferred_state[I]
            constants.reactive_indices[index] = LinearIndices(sz)[CartesianIndex.(I,J)]
        end
    else
        for (index, uniq) in enumerate(uniq_compounds)
            constants.reactive_indices[index] = findall(constants.compoundNames .== uniq)[1]
        end
    end

    # Constants (Yield, eDonor, mu_max, maint)
    yield_mu_file = file["Yield-Mu"]
    bac_names = collect(skipmissing(yield_mu_file[2:end,1]))
    names_Ymu = collect(skipmissing(yield_mu_file[1,2:end-1]))
    values_Ymu = reshape(collect(skipmissing(yield_mu_file[2:end,2:end-1])), length(bac_names), length(names_Ymu))

    if bac_names != constants.speciesNames
        throw(ErrorException("Bacteria species do not have the same name or are not in the same order"))
    end
    
    yield = values_Ymu[:, findall(names_Ymu .== "Yield")]
    eD = values_Ymu[:, findall(names_Ymu .== "eD")]

    constants.maintenance = values_Ymu[:,findall(names_Ymu .== "Maintenance")]
    constants.mu_max = values_Ymu[:,findall(names_Ymu .== "Max growth rate")]

    # Test is any of the maintenance or mu is unknown (missing)
    if ismissing(any(constants.maintenance .== 1)) || ismissing(any(constants.mu_max .== 1)) # The 1 is arbitrary, but needed a comparison to yield a "missing"
        constants.maintenance = NaN
        constants.mu_max = NaN
        println("Maintenance and maximum growth rate are not set, thus calculating dynamically. \nPlease make sure the equations and species match up in the code.\n")
    end

    # Constants (ReactionMatrix)
    bac_names_RM = collect(skipmissing(file["ReactionMatrix"][1, 2:end]))
    reaction_names = collect(skipmissing(file["ReactionMatrix"][2,2:end]))
    compounds = collect(skipmissing(file["ReactionMatrix"][3:end,1]))
    values_reacM = reshape(collect(skipmissing(file["ReactionMatrix"][3:end, 2:end])), length(compounds), length(reaction_names))

    first_compoundindex = findall(constants.compoundNames[1] .== compounds)
    last_compoundindex = findall(constants.compoundNames[end] .== compounds)

    if bac_names_RM != constants.speciesNames
        throw(ErrorException("Bacteria species do not have the same name or are not in the same order"))
    end


    if compounds[first_compoundindex[1]:last_compoundindex[1]] != constants.compoundNames
        throw(ErrorException("Compounds do not have the same name, or are not in the same order"))
    end

    iCat = findall(reaction_names .== "Cat")
    iAnab = findall(reaction_names .== "Anab")
    iDecay = findall(reaction_names .== "Decay")

    constants.MatrixMet = zeros(nCompounds+2, length(constants.speciesNames)) # +2 to include H2O and H
    constants.MatrixDecay = zeros(nCompounds+2, length(constants.speciesNames)) # +2 to include H2O and H

    for species = 1:length(constants.speciesNames)
        # get anabolism and catabolism
        cata = values_reacM[:, iCat[species]]
        ana = values_reacM[:, iAnab[species]]

        # get eDonor and Yield
        eD_species = eD[species]
        Y = yield[species]

        # Check whether eD has value -1 in cat
        eD_index = findall(constants.compoundNames .== eD_species)
        if cata[eD_index][1] != -1
            throw(ErrorException("Catabolism is not normalised towards the electron Donor"))
        end

        # Calculate metabolism matrix entry TODO: make sure this is correct!
        # It seems to me that if eD is also used in anabolism that this equation does not hold anymore!
        constants.MatrixMet[:, species] = cata ./ Y .+ ana

        # Set decay matrix entry
        constants.MatrixDecay[:, species] = values_reacM[:, iDecay[species]]
    end

    # Initialisation of bacteria
    settings.model_type = values_bac[names_bac .== "Initialisation method"][1]

    if settings.model_type in ("granule", "mature granule")
        bac_init.granule_radius = values_bac[names_bac .== "Starting granule radius"][1]
        bac_init.start_nBac = values_bac[names_bac .== "Starting number of bacteria (granule)"][1]

    elseif settings.model_type in ("suspension")
        bac_init.start_nColonies = values_bac[names_bac .== "Starting number of microcolonies (suspension)"][1]
        bac_init.start_nBacPerColony = values_bac[names_bac .== "Starting number of bacteria per microcolony (suspension)"][1]
    else
        throw(ErrorException("Initialisation method <$(settings.model_type)> is not a valid method."))
    end

    return [grid, bac_init, constants, settings, init_params]
end


# file_loc = string(Base.source_dir(), "\\","test_file.xlsx")
# grid, bac_init, constants, settings, init_params = loadPresetFile(file_loc);

# a = reshape([11;2;13;14;15], 5, :)
# b = reshape([10;20;12;4;20], :, 5)

# c = a .- b

# d = findmax(minimum(c, dims = 1))
a = range(1,10,10)
repeat(a, 1,2)

