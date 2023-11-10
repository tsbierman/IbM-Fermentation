function set_concentrations(conc_old, set_concs, mask)
    """
    Sets the initial concentration in the bioaggregate

    Conc is a matrix of nx*ny * compound for the concentration of each compound per grid cell
    set_concs is a vecstor of initial concentrations per compound in the bioaggregate 
    (nCompounds * 1)
    mask is a logical matrix (nx * ny) with cells for which the concentration needs to be set
    """

    set_concs = reshape(set_concs, 1, 1, :)
    conc = .!mask .* conc_old + mask .* set_concs
    
    return conc
end
