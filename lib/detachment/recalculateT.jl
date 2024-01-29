function recalculateT(T, x_index, y_index, kDet, grid, Visited, x_centre, y_centre)
    """
    This function recalculates T value at the gridcell[y_index, x_index] using 
    quadratic approximation of the gradient.

    Arguments
    T:                  A (ny,nx) matrix with the current detachment times
    y_index, x_index:   The coordinates of the gridcell that needs recalculation of the T value
    kDet:               The detachment constant determining how fast detachment takes place
    grid:               A "General" struct containing all parameters related to the grid
    Visited:            A (ny,nx) BitMatrix with per gridcell whether the correct T-Value has been calculated
    x_centre, y_centre: The coordinats of the centre of the granule

    Returns
    T_new:              A (ny, nx) matrix with the updated detachment times
    """

    # include(string(pwd(), "\\lib\\detachment\\calculateLocalDetachmentRate.jl"))
    # include(string(pwd(), "\\lib\\detachment\\computeRoot.jl"))

    right_visited = Visited[y_index, x_index + 1]
    left_visited = Visited[y_index, x_index - 1]

    if left_visited
        if right_visited                                                        # If both neighbouring T's are already Visited
            Tx = minimum([T[y_index, x_index - 1], T[y_index, x_index + 1]])    # take minimum of both
        else                                                                    # If only left is Visited
            Tx = T[y_index, x_index - 1]                                        # Copy its value
        end
    elseif right_visited                                                        # If only right is Visited
        Tx = T[y_index, x_index + 1]                                            # Copy its value
    else                                                                        # If neither are Visited
        Tx = Inf                                                                # Inf
    end

    top_visited = Visited[y_index + 1, x_index]
    bottom_visited = (y_index != 1) && Visited[y_index - 1, x_index]            # Why only check in this case

    if top_visited
        if bottom_visited                                                       # If both neighbouring T's are already Visited
            Ty = minimum([T[y_index - 1, x_index], T[y_index + 1, x_index]])    # take minimum of both
        else                                                                    # If only top is Visited
            Ty = T[y_index + 1, x_index]                                        # Copy its value
        end
    elseif bottom_visited                                                       # If only bottom is Visited
        Ty = T[y_index - 1, x_index]                                            # Copy its value
    else                                                                        # If neither are Visited
                                                                                # Inf
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
