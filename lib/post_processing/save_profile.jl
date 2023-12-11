function init_save_profile(constants, grid)
    """
    Initialise the results struct for the number of saves that are going
    to happen. Utilizes datatypes within the required precision with the lowest storage requirements.
    
    constants: struct containing all simulation constants
    grid: struct containing all information regarding the grid
    
    Results: structs containing all variables to be saved, initialised for the number of saves that are going to be made
    """
    include(string(pwd(), "\\lib\\Struct_Module"))

    nSaves = ceil(constants.simulation_end / constants.dT_save) + 1

    # Bacterial variables
    bac_saved = General()
    bac_saved.nBacs = zeros(UInt32, nSaves)                             # Vector, Unsigned Integer, 32 bit
    bac_saved.x = zeros(Float32, nSaves, constants.max_nBac)            # Matrix, Float 32 bit (single precision)
    bac_saved.y = zeros(Float32, nSaves, constants.max_nBac)            # Matrix, Float 32 bit (single precision)
    bac_saved.radius = zeros(Float32, nSaves, constants.max_nBac)       # Matrix, Float 32 bit (single precision)
    bac_saved.species = zeros(UInt8, nSaves, constants.max_nBac)        # Matrix, Unsigned Integer, 8 bit
    bac_saved.active = zeros(Bool, nSaves, constants.max_nBac)          # Matrix, Boolean
    bac_saved.mu = zeros(Float32, nSaves, constants.max_nBac)           # Matrix, Float 32 bit (single precision)

    # Concentration variable
    nCompounds = length(constants.compoundNames)
    conc_saved = zeros(Float32, nSaves, grid.nx, grid.ny, nCompounds)            # Matrix, Float 32 bit (single precision)

    # pH variable
    pH_saved = zeros(Float32, nSaves, grid.nx, grid.ny)

    # reactor properties
    reactor_saved = General()
    reactor_saved.bulk_concs = zeros(Float32, nSaves, nCompounds)       # Matrix, Float 32 bit (single precision)
    reactor_saved.HRT = zeros(Float32, nSaves)                          # Vector, Float 32 bit (single precision)
    reactor_saved.granule_density = zeros(Float32, nSaves)              # Vector, Float 32 bit (single precision)

    return bac_saved, conc_saved, pH_saved, reactor_saved
end


function save_profile(bac, conc, bulk_concentrations, pH, invHRT, Time, grid, constants, directory)
    """
    Save important variables along the central axis of the bio-aggregate
    
    bac:                    struct containing all information regarding the bacteria
    conc:                   matrix containing all concentrations per grid cell as of (ix, iy, compound)
    bulk_concentrations:    vector of the bulk liquid concentration of all compounds
    pH:                     matrix containing the pH value per grid cell as (ix, iy)
    Time:                   simulation time
    grid:                   struct containing all information regarding the grid
    constants:              struct containing all simulation constants
    directory:              directory where results are to be stored in
    """

    # Initialise of load previous values
    results_file = string(directory, "\\results2D.jdl2")
    if Time == 0
        bac_saved, conc_saved, pH_saved, reactor_saved = init_save_profile(constants, grid)
    else
        bac_saved, conc_saved, pH_saved, reactor_saved = load(results_file, "bac_saved", "conc_saved", "pH_saved", "reactor_saved")
    end
    
    # Set values
    iSave = ceil(Time  / constants.dT_save)  + 1

    # Bacterial variables
    nBacs = length(bac.x)
    bac_saved.nBacs[iSave] = nBacs
    bac_saved.x[iSave, 1:nBacs] = bac.x
    bac_saved.y[iSave, 1:nBacs] = bac.y
    bac_saved.radius[iSave, 1:nBacs] = bac.radius
    bac_saved.species[iSave, 1:nBacs] = bac.species
    bac_saved.active[iSave, 1:nBacs] = bac.active
    bac_saved.mu[iSave, 1:nBacs] = bac.mu

    # concentration variable
    conc_saved[iSave, :, :, :] = conc   # Save all concentrations

    # pH variable
    pH_saved[iSave, :, :] = pH          # Save all pH values

    # reactor_properties
    reactor_saved.bulk_concs[iSave, :] = bulk_concentrations
    reactor_saved.HRT[iSave] = 1 / invHRT
    reactor_saved.granule_density[iSave] = sum(bac.molarMass .* constants.bac_MW) / ((maximum(bac.y) - minimum(bac.y)) * (maximum(bac.x) - minimum(bac.x)) * 1e-6) # [g/m3], why use radius of bacteria and not diameter?

    # Save structs to file
    save(results_file, "bac_saved", bac_saved, "conc_saved", conc_saved, "pH_saved", pH_saved, "reactor_saved", reactor_saved)

end
