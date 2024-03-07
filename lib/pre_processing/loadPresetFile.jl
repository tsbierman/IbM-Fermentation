function loadPresetFile(filename)
    """
    This function reads the excel and stores the parameters in structs

    Arguments
    filename:           The excel file that will be read

    Returns
    grid:               A "General" struct containing all parameters related to the grid
    bac_init:           A "General" struct containing parameters need for initialisation of bacteria
    constants:          A "General" struct containing all the simulation constants
    settings:           A "General" struct containing all the settings of the simulation
    init_params:        A "General" struct containing the parameters values at the start of the simulation
    """
    # Initialisation of the structs
    grid_float = Float_struct()
    grid_int = Int_struct()
    constants_float = Float_struct()
    constants_vecfloat = VectorFloat_struct()
    constants_vecint = VectorInt_struct()
    constants_vecstring = VectorString_struct()
    constants_vecbool = VectorBool_struct()
    constants_matfloat = MatrixFloat_struct()
    settings_bool = Bool_struct()
    settings_string = String_struct()
    init_params = VectorFloat_struct()
    bac_init_float = Float_struct()
    bac_init_int = Int_struct()

    #Read file
    file = XLSX.readxlsx(filename)

    # initialise grid
    names_discr, values_discr = collect(skipmissing(file["Discretization"][:,1])), collect(skipmissing(file["Discretization"][:,2]))
    grid_float.dx = values_discr[names_discr .== "dx"][1]                                         # [m]
    grid_float.dy = values_discr[names_discr .== "dy"][1]                                         # [m]
    grid_float.dz = values_discr[names_discr .== "dz"][1]                                         # [m]
    grid_int.nx = values_discr[names_discr .== "nx"][1]                                         # [-]
    grid_int.ny = values_discr[names_discr .== "ny"][1]                                         # [-]
    grid_float.blayer_thickness = values_discr[names_discr .== "Boundary layer thickness"][1]     # [m]

    constants_float.Vg = (grid_float.dx ^ 3) * 1000                                                     # [L] Conversion m3 --> L
    constants_float.max_granule_radius = ((grid_int.nx - 4) * grid_float.dx) / 2                            # [m]

    # initialise constants (Time)
    constants_float.simulation_end = values_discr[names_discr .== "Simulation end"][1]            # [h]
    constants_float.dT = values_discr[names_discr .== "Initial dT diffusion"][1]                  # [h]
    constants_float.dT_bac = values_discr[names_discr .== "Initial dT bacteria"][1]               # [h]
    constants_float.dT_save = values_discr[names_discr .== "dT save"][1]                          # [h]
    constants_float.dT_analyse = values_discr[names_discr .== "dT analyse"][1]                    # [h]
    constants_float.dT_backup = values_discr[names_discr .== "dT backup"][1]                      # [h]

    settings_bool.dynamicDT = values_discr[names_discr .== "Dynamic dT"][1]                      # [Bool]

    # Only if dynamic time stepping is enabled
    if settings_bool.dynamicDT
        constants_vecint.nIterThreshold = values_discr[names_discr .== "nIterThreshold"]                            # [-]
        constants_vecint.iterThresholdDecrease = values_discr[names_discr .== "iterThresholdDecrease"]              # [-]
        constants_vecint.iterThresholdIncrease = values_discr[names_discr .== "iterThresholdIncrease"]              # [-]
        constants_float.initRESThresholdIncrease = values_discr[names_discr .== "initial RES threshold increase"][1]  # [mol/L/h]
        constants_vecint.nItersCycle = values_discr[names_discr .== "nIters per cycle"]                             # [-]
        constants_float.tolerance_no_convergence = values_discr[names_discr .== "tolerance no-convergence"][1]        # [mol/L/h]
        constants_float.maxRelDiffBulkConc = values_discr[names_discr .== "maximum relative bulk conc change"][1]     # [-]
        
        constants_float.maxDT = values_discr[names_discr .== "Maximum dT diffusion"][1]          # [h]
        constants_float.minDT = values_discr[names_discr .== "Minimum dT diffusion"][1]          # [h]
        constants_float.maxDT_bac = values_discr[names_discr .== "Maximum dT bacteria"][1]       # [h]
        constants_float.minDT_bac = values_discr[names_discr .== "Minimum dT bacteria"][1]       # [h]
    end

    # Constants (Diffusion)
    names_diff, values_diff = collect(skipmissing(file["Diffusion"][:,1])), collect(skipmissing(file["Diffusion"][:,2]))
    constants_vecstring.compoundNames = names_diff
    nCompounds = length(constants_vecstring.compoundNames)
    constants_vecfloat.diffusion_rates = values_diff                                                 # [m2/h]

    # Constants (Operational parameters)
    names_para, values_para = collect(skipmissing(file["Parameters"][:,1])), collect(skipmissing(file["Parameters"][:,2])) # It takes some extra empty rows, this removes that
    constants_float.pHsetpoint = values_para[names_para .== "pH setpoint"][1]                     # [-]
    constants_float.T = values_para[names_para .== "Temperature (K)"][1]                          # [K]
    constants_float.Vr = values_para[names_para .== "Representative volume"][1] * 1000            # [L]
    constants_float.reactor_density = values_para[names_para .== "Density reactor"][1]            # [g/L]

    settings_bool.variableHRT = values_para[names_para .== "Variable HRT"][1]                    # [Bool]
    init_params.invHRT = 1 ./ values_para[names_para .== "HRT"]                              # [1/h]

    # Only if varaible HRT is enabled
    if settings_bool.variableHRT
        constants_float.bulk_setpoint = values_para[names_para .== "Setpoint"][1]              # [mol/L]
        compound_name = values_para[names_para .== "Compound setpoint"][1] 
        constants_vecint.setpoint_index = findall(constants_vecstring.compoundNames .== compound_name) # [Vec(Int)]
    end

    # Constants (Bacteria)
    names_bac, values_bac = collect(skipmissing(file["Bacteria"][:,1])), collect(skipmissing(file["Bacteria"][:,2]))  # It takes some extra empty rows, this removes that
    constants_float.bac_MW = values_bac[names_bac .== "Molecular weight bacterium"][1]            # [g/mol]
    constants_float.bac_rho = values_bac[names_bac .== "Density bacterium"][1]                    # [g/m3]
    constants_vecint.max_nBac = values_bac[names_bac .== "Maximum nBacteria"]                   # [-] Maybe find better way of calculating this
    constants_vecbool.inactivationEnabled = values_bac[names_bac .== "Inactivation enabled"]     # [Bool]
    constants_float.min_bac_mass_grams = values_bac[names_bac .== "Minimum mass bacterium"][1]    # [g]
    constants_float.max_bac_mass_grams = values_bac[names_bac .== "Maximum mass bacterium"][1]    # [g]
    constants_float.bac_max_radius = values_bac[names_bac .== "Maximum radius bacterium"][1]      # [m]
    constants_float.kDist = values_bac[names_bac .== "kDist"][1]                                  # [-]
    constants_float.max_granule_radius = min(constants_float.max_granule_radius, values_bac[names_bac .== "Maximum granule radius"][1]) # [m] Either based on grid size or set value
    constants_float.kDet = values_bac[names_bac .== "Detachment constant"][1]                     # [1/m2.h]
    settings_string.detachment = values_bac[names_bac .== "Detachment method"][1]                  # [mechanistic, naive, none]


    # Constants (Solver)
    names_solv, values_solv = collect(skipmissing(file["Solver"][:,1])), collect(skipmissing(file["Solver"][:,2]))
    constants_float.diffusion_accuracy = values_solv[names_solv .== "Diffusion tolerance"][1]                 # [-]
    constants_float.steadystate_tolerance = values_solv[names_solv .== "Steady state RES threshold"][1]       # [mol/L/h]
    tol_abs = values_solv[names_solv .== "Concentration tolerance"][1]                                  # [mol/L]
    constants_float.correction_concentration_steady_state = tol_abs / constants_float.steadystate_tolerance         # Not really used anymore
    constants_vecstring.RESmethod = values_solv[names_solv .== "RES determination method"]                     # [max, mean, norm]

    constants_vecint.nDiffusion_per_SScheck = values_solv[names_solv .== "nIters diffusion per SS check"]   # [-]

    settings_bool.pHbulkCorrection = values_solv[names_solv .== "pH bulk concentration corrected"][1]        # Boolean
    settings_bool.pHincluded = values_solv[names_solv .== "pH solving included"][1]                          # Boolean
    # Speciation is always included if pH is included, else, read Boolean from excel
    settings_bool.speciation = values_solv[names_solv .== "Speciation included"][1] || settings_bool.pHincluded   # Boolean

    if settings_bool.pHincluded
        constants_float.pHtolerance = values_solv[names_solv .== "pH solver tolerance"][1]                    # [-]
    else
        constants_float.pHtolerance = NaN
    end

    settings_bool.structure_model = values_solv[names_solv .== "Structure model"][1]                         # Boolean
    if settings_bool.structure_model
        settings_string.type = values_solv[names_solv .== "Structure model type"][1]                           # [Neut, Comm, Comp, Copr]
    end

    settings_bool.parallelized = values_solv[names_solv .== "Parallelisation"][1]                            # Boolean

    # Constants (influent)
    names_temp, values_temp, condition_type_temp = file["Influent"][:,1], file["Influent"][:,2], file["Influent"][:,4] 
    names_infl, values_infl, condition_type = collect(skipmissing(names_temp)), collect(skipmissing(values_temp)), collect(skipmissing(condition_type_temp)) # It takes some extra empty rows, this removes that

    if names_infl != constants_vecstring.compoundNames
        throw(ErrorException("Compounds do not have the same name or are not in the same order"))
    end

    constants_vecbool.Dir_k = condition_type .== "D"
    constants_vecfloat.influent_concentrations = values_infl                                  # [mol/L]

    # Constants (initial conditions)
    names_init, values_init = collect(skipmissing(file["Initial condition"][:,1])), collect(skipmissing(file["Initial condition"][:,2]))

    if names_init != constants_vecstring.compoundNames
        throw(ErrorException("Compounds do not have the same name or are not in the same order"))
    end

    # For the species with a Dirichlet, the influent concentration is taken as the bulk concentration
    init_params.init_concs = copy(values_init)                                            # [mol/L]
    init_params.init_bulk_conc = copy(values_init)                                        # [mol/L]
    init_params.init_bulk_conc[constants_vecbool.Dir_k] = constants_vecfloat.influent_concentrations[constants_vecbool.Dir_k]

    # Constants (Equilibrium constants & charge matrix)
    thermodynamic_parameters = file["ThermoParam"]
    names_thermo = thermodynamic_parameters[:,1]
    nColumns = length(file["ThermoParam"][2,:]) - 4             # -4 due to Compound names, preferred subspecies, compound phase and cell with text, bit Hardcoded but should remain the same amount of columns

    # Divide in sections and locate H2O and H indices
    section_starts = findall(names_thermo .== constants_vecstring.compoundNames[1])
    section_ends = findall(names_thermo .== constants_vecstring.compoundNames[end])
    H2O_index = findall(names_thermo .== "H2O")
    H_index = findall(names_thermo .== "H")
    dG_rows = collect(section_starts[1]:section_ends[1])

    # Check the order of the compounds names
    compound_names_matrix = reshape(constants_vecstring.compoundNames, length(constants_vecstring.compoundNames),:) # Change shape to compare with thermodynamic names
    if thermodynamic_parameters[:,1][dG_rows] != (compound_names_matrix)
        throw(ErrorException("Compounds do not have the same name or are not in the same order"))
    end

    # Combine a section with the H2O and H indices
    dG_rows = push!(vec(dG_rows), H2O_index[1], H_index[1])

    preferred_state = thermodynamic_parameters[:,7][dG_rows]

    # Get next section indices
    Keq_start = section_starts[2][1]
    Keq_end = section_ends[2][1] + 2    # For adding H20 and H    
    constants_matfloat.Keq = thermodynamic_parameters[Keq_start:Keq_end, 2:5]

    # And the last section indices
    charge_start = section_starts[3][1]
    charge_end = section_ends[3][1] 
    charge_rows = collect(section_starts[3]:section_ends[3])

    # Check order of compounds
    if thermodynamic_parameters[:,1][charge_rows] != compound_names_matrix
        throw(ErrorException("Compounds do not have the same name or are not in the same order"))
    end

    # Extract charge matrix and substitute everything without charge
    charge_end = charge_end + 2     # For adding H2O and H
    temp_chrM = thermodynamic_parameters[charge_start:charge_end, 2:6]
    idx = findall(temp_chrM .== "NA")
    temp_chrM[idx] .= 0
    constants_matfloat.chrM = temp_chrM

    # Constants (Ks and Ki)
    ks_file = file["Ks"]
    ki_file = file["Ki"]

    # Extract and check microbial species Names
    constants_vecstring.speciesNames = collect(skipmissing(ks_file[:,1][2:end,1]))

    if collect(skipmissing(ki_file[:,1][2:end,1])) != constants_vecstring.speciesNames
        throw(ErrorException("Bacteria species do not have the same name or are not in the same order"))
    end

    # Get names of compounds, make a unique set and extract matrix of values
    comp_names_ks = collect(skipmissing(ks_file[1,1:end-1]))
    comp_names_ki = collect(skipmissing(ki_file[1,1:end-1]))
    uniq_compounds = unique([comp_names_ks; comp_names_ki])
    values_ks = reshape(collect(skipmissing(ks_file[2:end,2:end-1])),length(constants_vecstring.speciesNames), length(comp_names_ks))
    values_ki = reshape(collect(skipmissing(ki_file[2:end,2:end-1])),length(constants_vecstring.speciesNames), length(comp_names_ki))

    constants_matfloat.Ks = zeros(length(constants_vecstring.speciesNames),length(uniq_compounds))
    constants_matfloat.Ki = zeros(length(constants_vecstring.speciesNames),length(uniq_compounds))

    # Loop to check whether a compound has a Ks and Ki, if yes, store value
    for (index, uniq) in enumerate(uniq_compounds)
        columnidx = findall(uniq .== comp_names_ks)
        if length(columnidx) != 0
            constants_matfloat.Ks[:,columnidx] = values_ks[:,columnidx]
        end

        columnidx = findall(uniq .== comp_names_ki)
        if length(columnidx) != 0
            constants_matfloat.Ki[:,columnidx] = values_ki[:,columnidx]
        end
    end

    # Reactive indices
    constants_vecint.reactive_indices = zeros(Int, length(uniq_compounds))
    if settings_bool.speciation
        sz = (nCompounds+2, nColumns)

        for (index, uniq) in enumerate(uniq_compounds)
            I = findall(constants_vecstring.compoundNames .== uniq)[1][1]
            J = preferred_state[I]
            constants_vecint.reactive_indices[index] = LinearIndices(sz)[CartesianIndex.(I,J)]
        end
    else
        for (index, uniq) in enumerate(uniq_compounds)
            constants_vecint.reactive_indices[index] = findall(constants_vecstring.compoundNames .== uniq)[1]
        end
    end

    # Constants (Yield, eDonor, mu_max, maint)
    yield_mu_file = file["Yield-Mu"]
    bac_names = collect(skipmissing(yield_mu_file[2:end,1]))
    names_Ymu = collect(skipmissing(yield_mu_file[1,2:end-1]))
    values_Ymu = reshape(collect(skipmissing(yield_mu_file[2:end,2:end-1])), length(bac_names), length(names_Ymu))

    if bac_names != constants_vecstring.speciesNames
        throw(ErrorException("Bacteria species do not have the same name or are not in the same order"))
    end
    
    yield = values_Ymu[:, findall(names_Ymu .== "Yield")]
    eD = values_Ymu[:, findall(names_Ymu .== "eD")]

    constants_vecfloat.maintenance = dropdims(values_Ymu[:,findall(names_Ymu .== "Maintenance")], dims=2)
    constants_vecfloat.mu_max = dropdims(values_Ymu[:,findall(names_Ymu .== "Max growth rate")], dims=2)

    # Test whether any of the maintenance or mu is unknown (missing)
    if any(ismissing.(constants_vecfloat.maintenance)) || any(ismissing.(constants_vecfloat.mu_max))
        constants_vecfloat.maintenance = [NaN]
        constants_vecfloat.mu_max = [NaN]
        println("Maintenance and maximum growth rate are not set, thus calculating dynamically. \nPlease make sure the equations and species match up in the code.\n")
    end

    # Constants (ReactionMatrix)
    bac_names_RM = collect(skipmissing(file["ReactionMatrix"][1, 2:end]))
    reaction_names = collect(skipmissing(file["ReactionMatrix"][2,2:end]))
    compounds = collect(skipmissing(file["ReactionMatrix"][3:end,1]))
    values_reacM = reshape(collect(skipmissing(file["ReactionMatrix"][3:end, 2:end])), length(compounds), length(reaction_names))

    first_compoundindex = findall(constants_vecstring.compoundNames[1] .== compounds)
    last_compoundindex = findall(constants_vecstring.compoundNames[end] .== compounds)

    if bac_names_RM != constants_vecstring.speciesNames
        throw(ErrorException("Bacteria species do not have the same name or are not in the same order"))
    end

    if compounds[first_compoundindex[1]:last_compoundindex[1]] != constants_vecstring.compoundNames
        throw(ErrorException("Compounds do not have the same name, or are not in the same order"))
    end

    iCat = findall(reaction_names .== "Cat")
    iAnab = findall(reaction_names .== "Anab")
    iDecay = findall(reaction_names .== "Decay")

    constants_matfloat.MatrixMet = zeros(nCompounds, length(constants_vecstring.speciesNames))
    constants_matfloat.MatrixDecay = zeros(nCompounds, length(constants_vecstring.speciesNames))

    for species = 1:length(constants_vecstring.speciesNames)
        # get anabolism and catabolism
        cata = values_reacM[first_compoundindex[1]:last_compoundindex[1], iCat[species]]
        ana = values_reacM[first_compoundindex[1]:last_compoundindex[1], iAnab[species]]

        # get eDonor and Yield
        eD_species = eD[species]
        Y = yield[species]

        # Check whether eD has value -1 in cat
        eD_index = findall(constants_vecstring.compoundNames .== eD_species)
        if cata[eD_index][1] != -1
            throw(ErrorException("Catabolism is not normalised towards the electron Donor"))
        end

        # Calculate metabolism matrix entry TODO: make sure this is correct!
        # It seems to me that if eD is also used in anabolism that this equation does not hold anymore!
        constants_matfloat.MatrixMet[:, species] = cata ./ Y .+ ana

        # Set decay matrix entry
        constants_matfloat.MatrixDecay[:, species] = values_reacM[first_compoundindex[1]:last_compoundindex[1], iDecay[species]]
    end

    # Initialisation of bacteria
    settings_string.model_type = values_bac[names_bac .== "Initialisation method"][1]

    if settings_string.model_type in ("granule", "mature granule")
        bac_init_float.granule_radius = values_bac[names_bac .== "Starting granule radius"][1]
        bac_init_int.start_nBac = values_bac[names_bac .== "Starting number of bacteria (granule)"][1]

    elseif settings_string.model_type in ("suspension",)
        bac_init_int.start_nColonies = values_bac[names_bac .== "Starting number of microcolonies (suspension)"][1]
        bac_init_int.start_nBacPerColony = values_bac[names_bac .== "Starting number of bacteria per microcolony (suspension)"][1]
    else
        throw(ErrorException("Initialisation method <$(settings_string.model_type)> is not a valid method."))
    end

    return grid_float, grid_int, bac_init_float, bac_init_int, constants_float, constants_vecfloat, constants_vecint, constants_vecstring, constants_vecbool, constants_matfloat, settings_bool, settings_string, init_params
end
