function recalculateT(T, x_index, y_index, kDet, grid, Visited, x_centre, y_centre)
    """
    Recalculate T value at the gridcell[y_index, x_index] using quadratic approximation of the gradient
    T is a matrix with the current detachment times
    y_index, x_index are the coordinates of the gridcell that needs recalculation of the T value
    kDet is the detachment constant determining how fast detachment takes place
    grid is a struct containing all information regarding the spacial discretization
    Visited is a logical matrix with per gridcell whetehr the correct T-Value has been calculated
    x_centre, y_centre are the centre of the granule
    """

    include(string(pwd(), "\\lib\\detachment\\calculateLocalDetachmentRate.jl"))
    include(string(pwd(), "\\lib\\detachment\\computeRoot.jl"))

    right_visited = Visited[y_index, x_index + 1]
    left_visited = Visited[y_index, x_index - 1]

    if left_visited
        if right_visited    # If both neighbouring T's are already definite
            Tx = mininum(T[y_index, x_index - 1], T[y_index, x_index + 1]) # take minimum of both
        else                # If only left is definite
            Tx = T[y_index, x_index - 1] # Copy its value
        end
    elseif right_visited # If only right is definite
        Tx = T[y_index, x_index + 1] # Copy its value
    else # If neither are definite
        Tx = Inf
    end

    top_visited = Visited[y_index + 1, x_index]
    bottom_visited = (y_index != 1) && Visited[y_index - 1, x_index] # Why only check in this case

    if top_visited
        if bottom_visited    # If both neighbouring T's are already definite
            Ty = mininum(T[y_index - 1, x_index], T[y_index + 1, x_index]) # take minimum of both
        else                # If only top is definite
            Ty = T[y_index + 1, x_index] # Copy its value
        end
    elseif bottom_visited # If only bottom is definite
        Ty = T[y_index - 1, x_index] # Copy its value
    else # If neither are definite
        Ty = Inf
    end

    if isinf(Tx) && isinf(Ty)
        error("All neighbours have infinite time of crossing")
    end

    Fdet = calculateLocalDetachmentRate(x_index, y_index, kDet, grid, x_centre, y_centre)

    if Fdet == 0
        warning("Detachment speed equals 0, thus infinite time of crossing") # Fdet can never be 0 due to formula
        T_new = Inf
        return T_new
    end

    T_new = computeRoot(Tx, Ty, Fdet, grid.dx)
    return T_new
end
