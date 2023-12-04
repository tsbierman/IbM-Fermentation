function smoothing(phi, rhs, L_lhs)
    """
    Apply a Jacobian smoothing iteration on the system using convolution
    This uses u(k+1) = (f - A*u(k)) / A[2,2]

    phi: the matrix of what has to be smoothed (u in the equation)
        can be either concentrations or the error that will be smoothed
    rhs: the right hand side of the equation (f in the equation above)
        can be either the rhs of the diffusion-reaction or residuals
    L_lhs: A stencil (kernel) that is part of the left-hand-side (A in the equation)
    L_lhs consists of the do-nothing kernel - a correction of the adjusted Laplacian kernel (L0)
    So, L_lhs is [I] - (1/2)^(2*depth) * alpha * L
    """
    L_sm = deepcopy(L_lhs)
    L_sm[2,2] = 0
    phi = (rhs .- conv(phi, L_sm)[2:end-1,2:end-1]) ./ L_lhs[2,2]
    return phi
end
