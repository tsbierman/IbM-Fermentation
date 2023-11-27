function create_dirichlet_boundary(x, value)
    """
    Constructs a dirichlet boundary around the given matrix x
    x: matrix (:,:) (for 1 compound) around which to construct the boundary
    value: boundary value for this compound

    x_dirichlet: padded matrix x with the boundary condition
    """

    x_dirichlet = zeros(size(x) .+ 2)
    x_dirichlet[2:end-1, 2:end-1] = x
    x_dirichlet[[1, end], :] .= value
    x_dirichlet[:, [1, end]] .= value

    return x_dirichlet
end
