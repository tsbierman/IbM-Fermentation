function restriction(res_fine, L_restriction)
    """
    This function converts a matrix to a coarser grid, based on full weighting restriction

    Arguments
    res_fine:           A (n, m) matrix with the residual values
    L_restriction:      The restriction stencil (kernel)

    Returns
    res_coarse:         A ((n+1)/2, (m+1)/2) matrix with the restriced residual values
    """
    
    res_coarse = conv(res_fine, L_restriction)[2:end-1,2:end-1]
    res_coarse = res_coarse[1:2:end, 1:2:end]
    return res_coarse
end
