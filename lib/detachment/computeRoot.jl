function computeRoot(Tx, Ty, Fdet, dx)
    """
    Solve the quadratic equation for the gradient approximation of detachment front
    
    Tx, Ty: Time of detachment value of neighbouring gridcel in the x and y directions
    Fdet: speed of detachment at the specific gridcell
    dx: discretization resolution

    Returns: solution of the quadratic equation.
    """

    a = 0
    b = 0
    c = -(dx/Fdet)^2

    # if Tx is finite, add terms to parameters
    if isfinite(Tx)
        a = a + 1
        b = b - 2*Tx
        c = c + Tx^2
    end

    # Same for Ty
    if isfinite(Ty)
        a = a + 1
        b = b - 2*Ty
        c = c + Ty^2
    end

    # If all Infinite, return inf
    if a == 0
        root = Inf
        return root
    end

    D = sqrt(b^2 - 4 * a * c)

    if D < 0
        error("ValueError: Should always be above 0...")
    end

    # Positive solution is only valid
    root = (-b+D)/ (2*a)
    return root
end
