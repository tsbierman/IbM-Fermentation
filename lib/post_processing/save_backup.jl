function save_backup(bac, bulk_concs, invHRT, conc, reaction_matrix, pH, directory)
    """
    Save important variables required for restart at this point in time
    Will overwrite the last backup in order to always have the latest file
    
    bac:                struct containing all infromation regarding bacteria
    conc:               matrix containing all concentrations per grid cell (ny,nx,icompound)
    bulk_concs:         vector of bulk liquid concentrations of all compounds
    pH:                 matrix containing the pH value per grid cell (ny, nx)
    invHRT:             current 1/HRT value [1/h]
    reaction_matrix:    matrix containing all reaction rates per grid cell and compound (ny, nx, compound) [mol/L/h]
    directory:          directory where results are to be stored in
    """

    # Name file
    results_file = string(directory, "\\backup.jdl2")

    # Save everything
    save(results_file, "bac", bac, "bulk_concs", bulk_concs, "invHRT", invHRT, "conc", conc, "reaction_matrix", reaction_matrix, "pH", pH)
end
