function save_backup(bac_vecfloat, bac_vecint, bac_vecbool, bulk_concs, invHRT, conc, reaction_matrix, pH, directory)
    """
    This function saves important variables required for a restart at this point in time
    It will overwrite the last backup in order to always have the latest file
    
    Arguments
    bac_XYZ:            A struct containing bacterial parameters
    bulk_concs:         A (ncompounds,) vector of the bulk liquid concentrations of all compounds
    invHRT:             The inverse of the HRT [1/h]
    conc:               A (ny, nx, ncompounds) matrix containing all concentrations per gridcell
    reaction_matrix:    A (ny, nx, ncompounds) matrix containing all reaction rates per gridcell and compound [mol/L/h]
    pH:                 A (ny, nx) matrix containing the pH value per grid cell
    directory:          The directory where results are to be stored in
    """

    # Name file
    results_file = string(directory, "\\backup.jld2")

    # Save everything
    save(results_file, "bac_vecfloat", bac_vecfloat, "bac_vecint", bac_vecint, "bac_vecbool", bac_vecbool, "bulk_concs", bulk_concs, "invHRT", invHRT, "conc", conc, "reaction_matrix", reaction_matrix, "pH", pH)
end
