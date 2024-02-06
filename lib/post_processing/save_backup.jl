function save_backup(bac, bulk_concs, invHRT, conc, reaction_matrix, pH, directory)
    """
    This function saves important variables required for restart at this point in time
    It will overwrite the last backup in order to always have the latest file
    
    Arguments
    bac:                A "General" struct containing all parameters related to the bacteria
    conc:               A (ny, nx, ncompounds) matrix containing all concentrations per gridcell
    bulk_concs:         A (ncompounds,) vector of the bulk liquid concentrations of all compounds
    pH:                 A (ny, nx) matrix containing the pH value per grid cell
    invHRT:             The current 1/HRT value [1/h]
    reaction_matrix:    A (ny, nx, ncompounds) matrix containing all reaction rates per gridcell and compound [mol/L/h]
    directory:          The directory where results are to be stored in
    """

    # Name file
    results_file = string(directory, "\\backup.jld2")

    # Save everything
    save(results_file, "bac", bac, "bulk_concs", bulk_concs, "invHRT", invHRT, "conc", conc, "reaction_matrix", reaction_matrix, "pH", pH)
end
