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
    names, values = file["Discretization"][:,1], file["Discretization"][:,2]
    grid.dx = values[names .== "dx"]                                        # [m]
    grid.dy = values[names .== "dy"]                                        # [m]
    grid.nx = values[names .== "nx"]                                        # [-]
    grid.ny = values[names .== "ny"]                                        # [-]
    grid.blayer_thickness = values[names .== "Boundary layer thickness"]    # [m]

    constants.Vg = (grid.dx ^ 3) * 1000                                     # [L]
    constants.max_granule_radius = ((grid.nx - 4) * grid.dx) / 2            # [m]

    # initialise constants (Time)
    constants.simulation_end = values[names .== "Simulation end"]           # [h]
    constants.dT = values[names .== "Initial dT diffusion"]                 # [h]
    constants.dT_bac = values[names .== "Initial dT bacteria"]              # [h]
    constants.dT_save = values[names .== "dT save"]                         # [h]
    constants.dT_backup = values[names .== "dT backup"]                     # [h]

    settings.dynamicDT = values[names .== "Dynamic dT"][1]                  # [Bool]

    # Only if dynamic time stepping is enabled
    if settings.dynamicdt
        dynamicDT = General()   # Extra structure to store in the other structure
        dynamicDT.nIterThreshold = values[names .== "nIterThreshold"]
        dynamicDT.iterThresholdDecrease = values[names .== "iterThresholdDecrease"]
        dynamicDT.iterThresholdIncrease = values[names .== "iterThresholdIncrease"]
        dynamicDT.initRESThresholdIncrease = values[names .== "initial RES threshold increase"]
        dynamicDT.nItersCycle = values[names .== "nIters per cycle"]
        dynamicDT.tolerance_no_convergence = values[names .== "tolerance no-convergence"]
        dynamicDT.maxRelDiffBulkConc = values[names .== "maximum relative bulk conc chnage"]
        
        dynamicDT.maxDT = values[names .== "Maximum dT diffusion"]          # [h]
        dynamicDT.minDT = values[names .== "Minimum dT diffusion"]          # [h]
        dynamicDT.maxDT_bac = values[names .== "Maximum dT bacteria"]       # [h]
        dynamicDT.minDT_bac = values[names .== "Minimum dT bacteria"]       # [h]

        constants.dynamicDT = dynamicDT
    end

    # Constants (Diffusion)
    names, values = file["Diffusion"][:,1], file["Diffusion"][:,2]
    constants.compoundNames = names
    nCompounds = length(constants.compoundNames)
    constants.diffusion_rates = values                                      # [m2/h]

    # Constants (Operational parameters)
    names_temp, values_temp = file["Parameters"][:,1], file["Parameters"][:,2]
    names, values = collect(skipmissing(names_temp)), collect(skipmissing(values_temp)) # It takes some extra empty rows, this removes that
    constants.pHsetpoint = values[names .== "pH setpoint"]                  # [-]
    constants.T = values[names .== "Temperature (K)"]                       # [K]
    constants.Vr = values[names .== "Representative volume"] * 1000         # [L]
    constants.reactor_density = values[names .== "Density reactor"]         # [g/L]

    settings.variableHRT = values[names .== "Variable HRT"][1]              # [Bool]
    init_params.invHRT = 1/values[names .== "HRT"]                          # [1/h]

    # Only if varaible HRT is enabled
    if settings.variableHRT
        constants.bulk_setpoint = values[names .== "Setpoint"]              # [mol/L]
        compound_name = values[names .== "Compound setpoint"][1] 
        constants.setpoint_index = findall(constants.compoundNames .== compound_name) # [CartasianIndex]
    end

    # Constants (Bacteria)
    names_temp, values_temp = file["Bacteria"][:,1], file["Bacteria"][:,2]  
    names, values = collect(skipmissing(names_temp)), collect(skipmissing(values_temp)) # It takes some extra empty rows, this removes that
    constants.bac_MW = values[names .== "Molecular weight bacterium"]       # [g/mol]
    constants.bac_rho = values[names .== "Density bacterium"]               # [g/m3]
    constants.max_nBac = values[names .== "Maximum nBacteria"]              # [-] Maybe find better way of calculating this
    constants.inactivationEnabled = values[names .== "Inactivation enabled"][1] # [Bool]
    constants.min_bac_mass_grams = values[names .== "Minimum mass bacterium"] # [g]
    constants.max_bac_mass_grams = values[names .== "Maximum mass bacterium"] # [g]
    constants.bac_max_radius = values[names .== "Maximum radius bacterium"] # [m]
    constants.kDist = values[names .== "kDist"]                             # [-]
    constants.max_granule_radius = min(constants.max_granule_radius, values[names .== "Maximum granule radius"]) # [m] Either based on grid size or set value
    constants.kDet = values[names .== "Detachment constant"]                # [1/m2.h]
    settings.detachment = values[names .== "Detachment method"][1]          # [mechanistic, naive, none]


    # Constants (Solver)
    names, values = file["Solver"][:,1], file["Solver"][:,2]
    constants.diffusion_accuracy = values[names .== "Diffusion tolerance"]              # [-]
    constants.steadystate_tolerance = values[names .== "Steady state RES threshold"]    # [mol/L/h]
    tol_abs = values[names .== "Concentration tolerance"]                               # [mol/L]
    constants.correction_concentration_steady_state = tol_abs/constants.steadystate_tolerance
    constants.RESmethod = values[names .== "RES determination method"][1]

    constants.nDiffusion_per_SScheck = values[names .== "nIters diffusion per SS check"]

    settings.pHbulkCorrection = values[names .== "pH bulk concentration corrected"][1]
    settings.pHincluded = values[names .== "pH solving included"][1]
    # Speciation is always included if pH is included, else, read Boolean from excel
    settings.speciation = values[names .== "Speciation included"][1] || settings.pHincluded

    if settings.pHincluded
        constants.pHtolerance = values[names .== "pH solver tolerance"]
    else
        constants.pHtolerance = NaN
    end

    settings.structure_model = values[names .== "Structure model"][1]
    if settings.structure_model
        settings.type = values[names .== "Structure model type"][1]
    end

    settings.parallelized = values[names .== "Parallelisation"][1]

    # Constants (influent)
    names_temp, values_temp, condition_type_temp = file["Influent"][:,1], file["Influent"][:,2], file["Influent"][:,4] 
    names, values = collect(skipmissing(names_temp)), collect(skipmissing(values_temp)), collect(skipmissing(condition_type_temp)) # It takes some extra empty rows, this removes that

    if names != constants.compoundNames
        throw(ErrorException("Compounds do not have the same name or are not in the same order"))
    end

    constants.Dir_k = condition_type .== "D"
    constants.influent_concentrations = values                                  # [mol/L]

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
        throw(ErrorException("Species do not have the same name or are not in the same order"))
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
    for (index, value) in enumerate(uniq_compounds)
        columnidx = findall(value .== comp_names_ks)
        if length(columnidx) != 0
            constants.Ks[:,columnidx] = values_ks[:,columnidx]
        end

        columnidx = findall(value .== comp_names_ki)
        if length(columnidx) != 0
            constants.Ki[:,columnidx] = values_ki[:,columnidx]
        end
    end

    return [grid, bac_init, constants, settings, init_params]
end



import XLSX # Not needed in this file
file_loc = string(Base.source_dir(), "\\","test_file.xlsx") # Not needed in this file
file = XLSX.readxlsx(file_loc)

file_ks = file["Ks"]
file_ki = file["Ki"]

speciesNames = collect(skipmissing(file_ks[:,1][2:end,1]))
comp_names_ks = collect(skipmissing(file_ks[1,1:end-1]))
comp_names_ki = collect(skipmissing(file_ki[1,1:end-1]))

x = unique([comp_names_ks;comp_names_ki])

ks_vals = reshape(collect(skipmissing(file_ks[2:end,2:end-1])),length(speciesNames), length(comp_names_ks))
ki_vals = reshape(collect(skipmissing(file_ki[2:end,2:end-1])),length(speciesNames), length(comp_names_ki))

temp = zeros(length(speciesNames),length(x))

for (index, value) in enumerate(x)
    columnidx = findall(value .== comp_names_ki)
    if length(columnidx) != 0
        temp[:,columnidx] = ki_vals[:,columnidx]
    end
end

temp

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
