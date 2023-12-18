function create_dirichlet_boundary(x, value)
    """
    This function constructs a dirichlet boundary around the given matrix x

    Arguements
    x:              Matrix (any,any) (for 1 compound) around which to construct the boundary
    value:          Boundary value for this compound

    Returns
    x_dirichlet:    Padded matrix (any+2,any+2) with the boundary value around matrix x
    """

    x_dirichlet = zeros(size(x) .+ 2)
    x_dirichlet[2:end-1, 2:end-1] = x
    x_dirichlet[[1, end], :] .= value
    x_dirichlet[:, [1, end]] .= value

    return x_dirichlet
end
