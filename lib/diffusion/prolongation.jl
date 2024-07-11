function prolongation(phi_coarse, L_prolongation, sz)
    """
    This function linearly interpolates the coarse values to the fine values using convolution

    Arguments
    phi_coarse:         A (n, m) matrix with phi-values (could also be errors depending on where in the V-cycle this is called)
    L_prolongation:     The prolongation stencil (kernel)
    sz:                 The size of the new phi_grid

    Returns
    phi_fine:           A (2n-1,2m-1) matrix with the interpolated phi-values
    """

    phi_fine = zeros(sz)
    phi_fine[1:2:end, 1:2:end] = phi_coarse
    phi_fine = imfilter(phi_fine, reflect(centered(L_prolongation)), Fill(0))

    return phi_fine
end
