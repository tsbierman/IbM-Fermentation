function init_save_profile(constants_float, constants_vecint, constants_vecstring, grid_int)
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

    nSaves = ceil(Int, constants_float.simulation_end / constants_float.dT_save) + 1
    max_nBac = constants_vecint.max_nBac[1]

    # Bacterial variables
    bac_saved_vecint = VectorInt_struct()
    bac_saved_matfloat = MatrixFloat_struct()
    bac_saved_matint = MatrixInt_struct()
    bac_saved_matbool = MatrixBool_struct()

    bac_saved_vecint.nBacs = zeros(Int32, nSaves)                             # Vector, Unsigned Integer, 32 bit
    bac_saved_matfloat.x = zeros(Float64, nSaves, max_nBac)            # Matrix, Float 32 bit (single precision)
    bac_saved_matfloat.y = zeros(Float64, nSaves, max_nBac)            # Matrix, Float 32 bit (single precision)
    bac_saved_matfloat.radius = zeros(Float64, nSaves, max_nBac)       # Matrix, Float 32 bit (single precision)
    bac_saved_matint.species = zeros(Int32, nSaves, max_nBac)        # Matrix, Unsigned Integer, 8 bit
    bac_saved_matbool.active = zeros(Bool, nSaves, max_nBac)          # Matrix, Boolean
    bac_saved_matfloat.mu = zeros(Float64, nSaves, max_nBac)           # Matrix, Float 32 bit (single precision)

    # Concentration variable
    nCompounds = length(constants_vecstring.compoundNames)
    conc_saved = zeros(Float32, nSaves, grid_int.nx, grid_int.ny, nCompounds)            # Matrix, Float 32 bit (single precision)

    # pH variable
    pH_saved = zeros(Float32, nSaves, grid_int.nx, grid_int.ny)

    # reactor properties
    reactor_saved_matfloat = MatrixFloat_struct()
    reactor_saved_vecfloat = VectorFloat_struct()
    reactor_saved_matfloat.bulk_concs = zeros(Float64, nSaves, nCompounds)       # Matrix, Float 32 bit (single precision)
    reactor_saved_vecfloat.HRT = zeros(Float64, nSaves)                          # Vector, Float 32 bit (single precision)
    reactor_saved_vecfloat.granule_density = zeros(Float64, nSaves)              # Vector, Float 32 bit (single precision)

    return bac_saved_vecint, bac_saved_matfloat, bac_saved_matint, bac_saved_matbool, conc_saved, pH_saved, reactor_saved_matfloat, reactor_saved_vecfloat
end


function save_profile(bac_vecfloat, bac_vecint, bac_vecbool, conc, bulk_concentrations, pH, invHRT, Time, grid_float, grid_int, constants_float, constants_vecint, constants_vecstring, directory)
    """
    This function saves important variables of the whole simulation domain
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

    # Initialise or load previous values
    results_file = string(directory, "\\results2D.jld2")
    if Time == 0
        bac_saved_vecint, bac_saved_matfloat, bac_saved_matint, bac_saved_matbool, conc_saved, pH_saved, reactor_saved_matfloat, reactor_saved_vecfloat = init_save_profile(constants_float, constants_vecint, constants_vecstring, grid_int)
    else
        bac_saved_vecint, bac_saved_matfloat, bac_saved_matint, bac_saved_matbool, conc_saved, pH_saved, reactor_saved_matfloat, reactor_saved_vecfloat = load(results_file, "bac_saved_vecint", "bac_saved_matfloat", "bac_saved_matint", "bac_saved_matbool", "conc_saved", "pH_saved", "reactor_saved_matfloat", "reactor_saved_vecfloat")
    end
    
    # Set values
    # iSave = ceil(Time  / constants_float.dT_save) + 1
    iSave = ceil(Int, (Time+0.01) / constants_float.dT_save)

    # Bacterial variables
    nBacs = length(bac_vecfloat.x)
    bac_saved_vecint.nBacs[iSave] = nBacs
    bac_saved_matfloat.x[iSave, 1:nBacs] = bac_vecfloat.x
    bac_saved_matfloat.y[iSave, 1:nBacs] = bac_vecfloat.y
    bac_saved_matfloat.radius[iSave, 1:nBacs] = bac_vecfloat.radius
    bac_saved_matint.species[iSave, 1:nBacs] = bac_vecint.species
    bac_saved_matbool.active[iSave, 1:nBacs] = bac_vecbool.active
    bac_saved_matfloat.mu[iSave, 1:nBacs] = bac_vecfloat.mu

    # concentration variable
    conc_saved[iSave, :, :, :] = conc   # Save all concentrations

    # pH variable
    pH_saved[iSave, :, :] = pH          # Save all pH values

    # reactor_properties
    reactor_saved_matfloat.bulk_concs[iSave, :] = bulk_concentrations
    reactor_saved_vecfloat.HRT[iSave] = 1 / invHRT
    reactor_saved_vecfloat.granule_density[iSave] = sum(bac_vecfloat.molarMass .* constants_float.bac_MW) ./ ((maximum(bac_vecfloat.y) .- minimum(bac_vecfloat.y)) .* (maximum(bac_vecfloat.x) .- minimum(bac_vecfloat.x)) .* grid_float.dz) # [g/m3]

    # Save structs to file
    save(results_file, "bac_saved_vecint", bac_saved_vecint, "bac_saved_matfloat", bac_saved_matfloat, "bac_saved_matint", bac_saved_matint, "bac_saved_matbool", bac_saved_matbool, "conc_saved", conc_saved, "pH_saved", pH_saved, "reactor_saved_matfloat", reactor_saved_matfloat, "reactor_saved_vecfloat", reactor_saved_vecfloat)

end
