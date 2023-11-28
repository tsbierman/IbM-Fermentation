function restriction(r, L_restriction)
    """
    Converts the matrix r to a coarser grid, based on full weighting restriction
    r: residual values on a fine-grid
    L_restriction: restriction stencil (kernel)

    for a size n of a dimension of r, the new grid will have a dimension of size (n+1)/2
    """
    rhs = conv(r, L_restriction)[2:end-1,2:end-1]
    rhs = rhs[1:2:end, 1:2:end]
    return rhs
end
