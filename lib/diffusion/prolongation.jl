function prolongation(phi_coarse, L_prolongation, sz)
    """
    This function linearly interpolates the coarse values to the fine values using convolution
    phi_coarse: phi-values on the coarse grid
    L_prolongation: prolongation stencil (kernel)
    sz: size of the new phi_grid

    for a size n of a dimension of phi_coarse, phi_fine will have a dimension with size n*2 - 1
    """

    phi_fine = zeros(sz)
    phi_fine[1:2:end, 1:2:end] = phi_coarse
    phi_fine = conv(phi_fine, L_prolongation)[2:end-1,2:end-1]

    return phi_fine
end
