function residual(phi, rhs, L_lhs)
    """
    This function calculates the residuals based on r = f - Av
    Here, Au = f with u the exact solution and v an approximation.

    Arguments
    phi:                A matrix with the phi-values (v in the equation)
    rhs:                The right-hand-side of the diffusion reaction equation (f in the equation)
    L_lhs:              A stencil (kernel) that is part of the left-hand-side (A in the equation)
                        L_lhs consists of the do-nothing kernel - a correction of the adjusted Laplacian kernel (L0)
                        So, L_lhs is [I] - (1/2)^(2*depth) * alpha * L

    Returns
    r:                  the calculated residuals
    """
    
    # r = rhs .- conv(phi, L_lhs)[2:end-1,2:end-1]
    r = rhs .- imfilter(phi, reflect(centered(L_lhs)), Fill(0))
    return r
end
