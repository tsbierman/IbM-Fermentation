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
    names_discr, values_discr = file["Discretization"][:,1], file["Discretization"][:,2]
    grid.dx = values_discr[names_discr .== "dx"]                                        # [m]
    grid.dy = values_discr[names_discr .== "dy"]                                        # [m]
    grid.nx = values_discr[names_discr .== "nx"]                                        # [-]
    grid.ny = values_discr[names_discr .== "ny"]                                        # [-]
    grid.blayer_thickness = values_discr[names_discr .== "Boundary layer thickness"]    # [m]

    constants.Vg = (grid.dx ^ 3) * 1000                                     # [L]
    constants.max_granule_radius = ((grid.nx - 4) * grid.dx) / 2            # [m]

    # initialise constants (Time)
    constants.simulation_end = values_discr[names_discr .== "Simulation end"]           # [h]
    constants.dT = values_discr[names_discr .== "Initial dT diffusion"]                 # [h]
    constants.dT_bac = values_discr[names_discr .== "Initial dT bacteria"]              # [h]
    constants.dT_save = values_discr[names_discr .== "dT save"]                         # [h]
    constants.dT_backup = values_discr[names_discr .== "dT backup"]                     # [h]

    settings.dynamicDT = values_discr[names_discr .== "Dynamic dT"][1]                  # [Bool]

    # Only if dynamic time stepping is enabled
    if settings.dynamicdt
        dynamicDT = General()   # Extra structure to store in the other structure
        dynamicDT.nIterThreshold = values_discr[names_discr .== "nIterThreshold"]
        dynamicDT.iterThresholdDecrease = values_discr[names_discr .== "iterThresholdDecrease"]
        dynamicDT.iterThresholdIncrease = values_discr[names_discr .== "iterThresholdIncrease"]
        dynamicDT.initRESThresholdIncrease = values_discr[names_discr .== "initial RES threshold increase"]
        dynamicDT.nItersCycle = values_discr[names_discr .== "nIters per cycle"]
        dynamicDT.tolerance_no_convergence = values_discr[names_discr .== "tolerance no-convergence"]
        dynamicDT.maxRelDiffBulkConc = values_discr[names_discr .== "maximum relative bulk conc chnage"]
        
        dynamicDT.maxDT = values_discr[names_discr .== "Maximum dT diffusion"]          # [h]
        dynamicDT.minDT = values_discr[names_discr .== "Minimum dT diffusion"]          # [h]
        dynamicDT.maxDT_bac = values_discr[names_discr .== "Maximum dT bacteria"]       # [h]
        dynamicDT.minDT_bac = values_discr[names_discr .== "Minimum dT bacteria"]       # [h]

        constants.dynamicDT = dynamicDT
    end

    # Constants (Diffusion)
    names_diff, values_diff = file["Diffusion"][:,1], file["Diffusion"][:,2]
    constants.compoundNames = names_diff
    nCompounds = length(constants.compoundNames)
    constants.diffusion_rates = values_diff                                      # [m2/h]

    # Constants (Operational parameters)
    names_temp, values_temp = file["Parameters"][:,1], file["Parameters"][:,2]
    names_para, values_para = collect(skipmissing(names_temp)), collect(skipmissing(values_temp)) # It takes some extra empty rows, this removes that
    constants.pHsetpoint = values_para[names_para .== "pH setpoint"]                  # [-]
    constants.T = values_para[names_para .== "Temperature (K)"]                       # [K]
    constants.Vr = values_para[names_para .== "Representative volume"] * 1000         # [L]
    constants.reactor_density = values_para[names_para .== "Density reactor"]         # [g/L]

    settings.variableHRT = values_para[names_para .== "Variable HRT"][1]              # [Bool]
    init_params.invHRT = 1/values_para[names_para .== "HRT"]                          # [1/h]

    # Only if varaible HRT is enabled
    if settings.variableHRT
        constants.bulk_setpoint = values_para[names_para .== "Setpoint"]              # [mol/L]
        compound_name = values_para[names_para .== "Compound setpoint"][1] 
        constants.setpoint_index = findall(constants.compoundNames .== compound_name) # [CartasianIndex]
    end

    # Constants (Bacteria)
    names_temp, values_temp = file["Bacteria"][:,1], file["Bacteria"][:,2]  
    names_bac, values_bac = collect(skipmissing(names_temp)), collect(skipmissing(values_temp)) # It takes some extra empty rows, this removes that
    constants.bac_MW = values_bac[names_bac .== "Molecular weight bacterium"]       # [g/mol]
    constants.bac_rho = values_bac[names_bac .== "Density bacterium"]               # [g/m3]
    constants.max_nBac = values_bac[names_bac .== "Maximum nBacteria"]              # [-] Maybe find better way of calculating this
    constants.inactivationEnabled = values_bac[names_bac .== "Inactivation enabled"][1] # [Bool]
    constants.min_bac_mass_grams = values_bac[names_bac .== "Minimum mass bacterium"] # [g]
    constants.max_bac_mass_grams = values_bac[names_bac .== "Maximum mass bacterium"] # [g]
    constants.bac_max_radius = values_bac[names_bac .== "Maximum radius bacterium"] # [m]
    constants.kDist = values_bac[names_bac .== "kDist"]                             # [-]
    constants.max_granule_radius = min(constants.max_granule_radius, values_bac[names_bac .== "Maximum granule radius"]) # [m] Either based on grid size or set value
    constants.kDet = values_bac[names_bac .== "Detachment constant"]                # [1/m2.h]
    settings.detachment = values_bac[names_bac .== "Detachment method"][1]          # [mechanistic, naive, none]


    # Constants (Solver)
    names_solv, values_solv = file["Solver"][:,1], file["Solver"][:,2]
    constants.diffusion_accuracy = values_solv[names_solv .== "Diffusion tolerance"]              # [-]
    constants.steadystate_tolerance = values_solv[names_solv .== "Steady state RES threshold"]    # [mol/L/h]
    tol_abs = values_solv[names_solv .== "Concentration tolerance"]                               # [mol/L]
    constants.correction_concentration_steady_state = tol_abs/constants.steadystate_tolerance
    constants.RESmethod = values_solv[names_solv .== "RES determination method"][1]

    constants.nDiffusion_per_SScheck = values_solv[names_solv .== "nIters diffusion per SS check"]

    settings.pHbulkCorrection = values_solv[names_solv .== "pH bulk concentration corrected"][1]
    settings.pHincluded = values_solv[names_solv .== "pH solving included"][1]
    # Speciation is always included if pH is included, else, read Boolean from excel
    settings.speciation = values_solv[names_solv .== "Speciation included"][1] || settings.pHincluded

    if settings.pHincluded
        constants.pHtolerance = values_solv[names_solv .== "pH solver tolerance"]
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

    # WAS HERE WITH CORRECTING NAMES



    # Constants (initial conditions)
    names, values = file["Initial condition"][:,1], file["Initial condition"][:,2]

    if names != constants.compoundNames
        throw(ErrorException("Compounds do not have the same name or are not in the same order"))
    end

    # For the species with a Dirichlet, the influent concentration is taken as the bulk concentration
    init_params.init_concs = values                                             # [mol/L]
    init_params.init_bulk_conc = init_params.init_concs                         # [mol/L]
    init_params.init_bulk_conc[constants.Dir_k] = constants.influent_concentrations[constants.Dir_k]

    # Constants (Equilibrium constants & charge matrix)
    thermodynamic_parameters = file["ThermoParam"]
    names = thermodynamic_parameters[:,1]
    nColumns = length(file["ThermoParam"][2,:])-4 # -4 due to Compound names, preferred subspecies, compound phase and cell with text 

    # Divide in sections and locate H2O and H indices
    section_starts = findall(names .== constants.compoundNames[1])
    section_ends = findall(names .== constants.compoundNames[end])
    H2O_index = findall(names .== "H2O")
    H_index = findall(names .== "H")
    dG_rows = collect(section_starts[1]:section_ends[1])

    # Check the order of the compounds names
    if thermodynamic_parameters[:,1][dG_rows] != constants.compoundNames
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
    if thermodynamic_parameters[:,1][charge_rows] != constants.compoundNames
        throw(ErrorException("Compounds do not have the same name or are not in the same order"))
    end

    # Extract charge matrix and substitute everything without charge
    carge_end = charge_end + 2 # For adding H2O and H
    constants.chrM = thermodynamic_parameters[Keq_start:Keq_end, 2:6]
    idx = findall(chrM .== "NA")
    chrM[idx] .= 0

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
        columnidx = findall(value .== comp_names_ks)
        if length(columnidx) != 0
            constants.Ks[:,columnidx] = values_ks[:,columnidx]
        end

        columnidx = findall(value .== comp_names_ki)
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
            constants.reactive_indices[index] = findall(constants.compoundNames .== uniq)
        end
    end

    # Constants (Yield, eDonor, mu_max, maint)
    yield_mu_file = file["Yield-Mu"]
    bac_names = yield_mu_file[2:end,1]
    names = yield_mu_file[1,2:end-1]
    values = yield_mu_file[2:end,2:end-1]

    if bac_names != constants.speciesNames
        throw(ErrorException("Bacteria species do not have the same name or are not in the same order"))
    end

    yield = values[:, findall(names .== "Yield")[1][2]]
    eD = values[:, findall(names .== "eD")[1][2]]

    constants.maintenance = values[:, findall(names .== "Maintenance")[1][2]]
    constants.mu_max = values[:, findall(names .== "Max growth rate")[1][2]]

    # Test is any of the maintenance or mu is unknown (missing)
    if ismissing(any(constants.maintenance .== 1)) || ismissing(any(constants.mu_max)) # The 1 is arbitrary, but needed a comparison to yield a "missing"
        constants.maintenance = NaN
        constants.mu_max = NaN
        println("Maintenance and maximum growth rate are not set, thus calculating dynamically. \nPlease make sure the equations and species match up in the code.\n")
    end

    # Constants (ReactionMatrix)
    bac_names = file["ReactionMatrix"][1, 2:end]
    bac_names = bac_names[1:3:end]
    reaction_names = file["ReactionMatrix"][2,2:end]
    values_reacM = file["ReactionMatrix"][3:end, 2:end]
    compounds = file["ReactionMatrix"][3:end,1]

    first_compoundindex = findall(constants.compoundNames[1] .== compounds)
    last_compoundindex = findall(constants.compoundNames[end] .== compounds)

    if bac_names != constants.speciesNames
        throw(ErrorException("Bacteria species do not have the same name or are not in the same order"))
    end

    if compounds[first_compoundindex:last_compoundindex] != constants.compoundNames
        throw(ErrorException("Compounds do not have the same name, or are not in the same order"))
    end

    iCat = findall(reaction_names .== "Cat")
    iAnab = findall(reaction_names .== "Anab")
    iDecay = findall(reaction_names .== "Decay")

    constants.MatrixMet = zeros(nCompounds+2, length(constants.speciesNames)) # +2 to include H2O and H
    constants.MatrixDecay = zeros(nCompounds+2, length(constants.speciesNames)) # +2 to include H2O and H

    for species = 1:length(constants.speciesNames)
        # get anabolism and catabolism
        cata = values_reacM[:, iCat[species][2]]
        ana = values_reacM[:, iAnab[species][2]]

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
        constants.MatrixMet[:, species] = cata/Y + ana

        # Set decay matrix entry
        constants.MatrixDecay[:, species] = values_reacM[:, iDecay[species][2]]
    end

    # Initialisation of bacteria
    names_temp, values_temp = file["Bacteria"][:,1], file["Bacteria"][:,2]  
    names, values = collect(skipmissing(names_temp)), collect(skipmissing(values_temp)) # It takes some extra empty rows, this removes that
    settings.model_type = values[names .== "Initialisation method"]

    if settings.model_type in ("granule", "mature granule")
        bac_init.granule_radius = values[names .== "Starting granule radius"]
        bac_init.start_nBac = values[names .== "Starting number of bacteria (granule)"]

    elseif settings.model_type in ("suspension")
        bac_init.start_nColonies = values[names .== "Starting number of microcolonies (suspension)"]
        bac_init.start_nBacPerColony = values[names .== "Starting number of bacteria per microcolony (suspension)"]
    else
        throw(ErrorException("Initialisation method <$(setting.model_type)> is not a valid method."))
    end

    return [grid, bac_init, constants, settings, init_params]
end



import XLSX # Not needed in this file
file_loc = string(Base.source_dir(), "\\","test_file.xlsx") # Not needed in this file
file = XLSX.readxlsx(file_loc)

bac_names = file["ReactionMatrix"][1, 2:end]
bac_names = bac_names[1:3:end]
values_reacM = file["ReactionMatrix"][3:end, 2:end]
reaction_names = file["ReactionMatrix"][2,2:end]

compounds = file["ReactionMatrix"][3:end,1]
cat_indices = findall(reaction_names .== "Cat")
ana_indices = findall(reaction_names .== "Anab")

cata = values_reacM[:, cat_indices[1][2]]
ana = values_reacM[:, ana_indices[1][2]]

names
# yield_mu_file = file["Yield-Mu"]
# names2 = yield_mu_file[1,2:end-1]
# values2 = yield_mu_file[2:end,2:end-1]
# x = values2[:,findall(names2 .== "eD")[1][2]]

# names2, values2 = file["Diffusion"]["A"], file["Diffusion"]["B"]
# compoundNames = names2

# eD_index = findall(compoundNames .== x[1])
# cata[eD_index][1] == -1


# for species = 1:length(compounds)
#     println(species)
# end



# mu = values[:,findall(names .== "Max growth rate")[1][2]]
# maint = values[:,findall(names .== "Maintenance")[1][2]]
# names2, values2 = file["Diffusion"]["A"], file["Diffusion"]["B"]
# compoundNames = names2
# length(compoundNames)

# preferred_state = file["ThermoParam"]["G2:G11"]

# I = findall(compoundNames .== "CO2")[1][1]
# J = preferred_state[I]
# sz = (10,5)
# index = LinearIndices(sz)[CartesianIndex.(I,J)]

# length(file["ThermoParam"][2,:])-4




# file_ks = file["Ks"]
# file_ki = file["Ki"]

# speciesNames = collect(skipmissing(file_ks[:,1][2:end,1]))
# comp_names_ks = collect(skipmissing(file_ks[1,1:end-1]))
# comp_names_ki = collect(skipmissing(file_ki[1,1:end-1]))

# x = unique([comp_names_ks;comp_names_ki])

# ks_vals = reshape(collect(skipmissing(file_ks[2:end,2:end-1])),length(speciesNames), length(comp_names_ks))
# ki_vals = reshape(collect(skipmissing(file_ki[2:end,2:end-1])),length(speciesNames), length(comp_names_ki))

# temp = zeros(length(speciesNames),length(x))

# for (index, value) in enumerate(x)
#     columnidx = findall(value .== comp_names_ki)
#     if length(columnidx) != 0
#         temp[:,columnidx] = ki_vals[:,columnidx]
#     end
# end

# temp

# names, values = file["ThermoParam"][:,8], file["Bacteria"][:,2]
# names
# # values[names .== "Simulation end"]

# thermodynamic_parameters= file["ThermoParam"]
# names2, values2 = file["Diffusion"]["A"], file["Diffusion"]["B"]
# compoundNames = names2
# names = thermodynamic_parameters["A"]
# section_starts = findall(names .== compoundNames[1])
# section_ends = findall(thermodynamic_parameters["A"] .== compoundNames[end])
# dG_rows = collect(section_starts[1][1]:section_ends[1][1])
# H2O_index = findall(names .== "H2O")
# H_index = findall(names .== "H")
# dG_rows = push!(vec(dG_rows), H2O_index[1][1], H_index[1][1])
# thermodynamic_parameters["H"][dG_rows]

# Keq_start = section_starts[3][1]
# Keq_end = section_ends[3][1] + 2
# Keq_rows = collect(section_starts[3][1]:section_ends[2][1])
# # Keq_rows = push!(vec(Keq_rows), H2O_index[2][1], H_index[2][1])
# chrM = thermodynamic_parameters[Keq_start:Keq_end, 2:6]
# idx = findall(chrM .== "NA")
# chrM[idx] .= 0
# chrM


# dG_rows = [dG_rows, H2O_index[1],H_index[1]]

# typeof(dG_rows[2])
# thermodynamic_parameters["G"][dG_rows]
# thermodynamic_parameters["A"][dG_rows] == compoundNames

# names, values, cond = file["Influent"]["A"], file["Influent"]["B"], file["Influent"]["D"]
# cond = collect(skipmissing(cond))
# cond

# names, values = file["Initial condition"]["A"], file["Initial condition"]["B"]
# names2, values2, condition_type = file["Influent"]["A1:A8"], file["Influent"]["B1:B8"], file["Influent"]["D1:D8"]
# Dir_k = condition_type .== "D"

# values[Dir_k] = values2[Dir_k]


# typeof(values[names .== "Compound setpoint"][1])


# names, values = file["Parameters"]["A1:A12"], file["Parameters"]["B1:B12"]
# names2, values2 = file["Diffusion"]["A"], file["Diffusion"]["B"]
# compoundNames = names2
# (compound_name = values[names .== "Compound setpoint"][1])
# (findall(compoundNames .== compound_name))
