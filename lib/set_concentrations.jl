function set_concentrations!(conc_old, set_concs, mask)
    """
    This function sets the concentrations, but only at the places indicated by mask

    Arguments
    conc:               A (ny, nx, nCompounds) matrix containing all concentrations per gridcell
    set_concs:          A (nCompounds,) vector of concentrations that need to be set
    mask:               A (ny,nx) BitMatrix indicating whether a gridcell should have the new concentrations or the old ones

    Returns
    conc                A (ny, nx, nCompounds) matrix containing all updated concentrations per gridcell
    """

    set_concs = reshape(set_concs, 1, 1, :)
    conc = .!mask .* conc_old + mask .* set_concs
    
    return conc
end
