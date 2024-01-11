function init_save_slice(constants, grid)
    """
    This function initialises the resulting structs that will contain saved information.
    It utilizes datatypes within the required precision with the lowest storage requirements.
    
    Arguments
    constants:              A "General" struct containing all the simulation constants
    grid:                   A "General" struct containing all parameters related to the grid
    
    Returns
    XXX_saved               A "General" struct that containing all variables of "General" struct XXX, that need to be saved
                            They are initialised for the number of saves that are going to be made
    """
    # include(string(pwd(), "\\lib\\Struct_Module"))

    nSaves = ceil(Int, constants.simulation_end / constants.dT_save) + 1

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
    conc_saved = zeros(Float32, nSaves, grid.nx, nCompounds)            # Matrix, Float 32 bit (single precision)

    # pH variable
    pH_saved = zeros(Float32, nSaves, grid.nx)

    # reactor properties
    reactor_saved = General()
    reactor_saved.bulk_concs = zeros(Float32, nSaves, nCompounds)       # Matrix, Float 32 bit (single precision)
    reactor_saved.HRT = zeros(Float32, nSaves)                          # Vector, Float 32 bit (single precision)
    reactor_saved.granule_density = zeros(Float32, nSaves)              # Vector, Float 32 bit (single precision)

    return bac_saved, conc_saved, pH_saved, reactor_saved
end


function save_slice(bac, conc, bulk_concentrations, pH, invHRT, Time, grid, constants, directory)
    """
    This function saves important variables along the central axis of the bio-aggregate
    It does so in other structs that are generated or loaded
    
    bac:                    A "General" struct containing all parameters related to the bacteria
    conc:                   A (ny, nx, ncompounds) matrix containing all concentrations per gridcell
    bulk_concentrations:    A (ncompounds,) vector of the bulk liquid concentrations of all compounds
    pH:                     A (ny, nx) matrix containing the pH value per grid cell
    Time:                   The simulation time (not the struct!)
    grid:                   A "General" struct containing all parameters related to the grid
    constants:              A "General" struct containing all the simulation constants
    directory:              The directory where results are to be stored in
    """
    # Initialise of load previous values
    results_file = string(directory, "\\results1D.jld2")
    if Time == 0
        bac_saved, conc_saved, pH_saved, reactor_saved = init_save_slice(constants, grid)
    else
        bac_saved, conc_saved, pH_saved, reactor_saved = load(results_file, "bac_saved", "conc_saved", "pH_saved", "reactor_saved")
    end
    
    # Set values
    iSave = ceil(Int, (Time+0.01) / constants.dT_save)

    if Time >= constants.simulation_end
        iSave = iSave + 1
    end

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
    conc_saved[iSave, :, :] = conc[ceil(Int, grid.ny/2), :, :]   # Save a horizontal slice through the centre of the granule (Only makes sense in granule mode)

    # pH variable
    pH_saved[iSave, :] = pH[ceil(Int, grid.ny/2), :]             # Save a horizontal slice through the centre of the granule (Only makes sense in granule mode)

    # reactor_properties
    reactor_saved.bulk_concs[iSave, :] = bulk_concentrations
    reactor_saved.HRT[iSave] = 1 / invHRT
    reactor_saved.granule_density[iSave] = sum(bac.molarMass .* constants.bac_MW) / ((maximum(bac.y) - minimum(bac.y)) * (maximum(bac.x) - minimum(bac.x)) * grid.dz) # [g/m3]

    # Save structs to file
    save(results_file, "bac_saved", bac_saved, "conc_saved", conc_saved, "pH_saved", pH_saved, "reactor_saved", reactor_saved)

end
