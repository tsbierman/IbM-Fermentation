function restriction(res_fine, L_restriction)
    """
    Converts the matrix r to a coarser grid, based on full weighting restriction
    res_fine: residual values on a fine-grid
    L_restriction: restriction stencil (kernel)

    for a size n of a dimension of r, the new grid will have a dimension of size (n+1)/2
    """
    res_coarse = conv(res_fine, L_restriction)[2:end-1,2:end-1]
    res_coarse = res_coarse[1:2:end, 1:2:end]
    return res_coarse
end
