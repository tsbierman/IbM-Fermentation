function getFreeNeighbourCount(i, j, Visited)
    """
    Count the number of non-biomass neighbouring gridcells for gridcel[j,i]
    """
    nFreeNb = 0
    if Visited[j, i+1] == 1
        nFreeNb = nFreeNb + 1
    end

    if Visited[j, i-1] == 1
        nFreeNb = nFreeNb + 1
    end

    if Visited[j+1, i] == 1
        nFreeNb = nFreeNb + 1
    end

    if Visited[j-1, i] == 1
        nFreeNb = nFreeNb + 1
    end
    return nFreeNb
end


function string_to_coordinates(string_coords)
    """
    Turns the string into an y and x- coordinate
    """
    comma_index = findfirst(",", string_coords)[1]
    y_coord = parse(Int64, string_coords[1:comma_index-1])
    x_coord = parse(Int64, string_coords[comma_index+1:end])
    return y_coord, x_coord
end


function calcTimeOfDetach(bac, grid, grid2bac, grid2nBacs, constants)
    """
    This function calculates the time of detachment for each gridcell in the simulation domain
    Bac is a struct that contains all information regarding the bacteria
    grid is a struct containing all information regarding spacial discretization
    grid2bac is a nx * nY * ? matrix that for each cell contains which bacteria is located 
    there. The number corresponds to the index in the bac struct
    grid2nBacs is a (nx * ny) matrix containing the number of bacteria in each grid cell

    It returns a matrix with per gridcell the time of detachment
    """
    include(string(pwd(), "\\lib\\determine_where_bacteria_in_grid.jl"))
    include(string(pwd(), "\\lib\\determine_diffusion_region.jl"))
    include(string(pwd(), "\\lib\\detachment\\calculateLocalDetachmentRate.jl"))
    include(string(pwd(), "\\lib\\detachment\\recalculateT.jl"))
    include(string(pwd(), "\\test\\plotLogicalGrid.jl")) # For testing visualisation

    # Extract variables
    kDet = constants.kDet

    # Calculate centre of granule
    x_centre = mean(bac.x)
    y_centre = mean(bac.y)

    # Where are bacteria located in the grid? The grid2nBacs is slightly extended and 
    # morphed into a logical matrix
    detachment_grid = deepcopy(grid) # To keep the original grid values from changing
    detachment_grid.blayer_thickness = constants.kDist * constants.bac_max_radius * 2 # Distance factor times diameter
    # Returns logical matrix with per cell whether it is in a certain region. As we set the the boundary layer very small,
    # this will only be the aggregate of the granule
    aggregate, _ = determine_diffusion_region(grid2bac, grid2nBacs, bac, detachment_grid)
    biofilm = aggregate .> 0
    # return plotLogicalGrid(detachment_grid, aggregate) # Plot test for the aggregate

    # Create matrices
    T = zeros(size(grid2nBacs))
    Visited = .!biofilm
    # return plotLogicalGrid(detachment_grid, Visited) # Plot test for Visited

    # Find narrow band
    kernel = zeros(3,3)
    kernel[[1,3], 2] .= -1/4
    kernel[2, [1,3]] .= -1/4
    kernel[2,2] = 1
    # Narrow_band is the outermost grids of the aggregate
    Narrow_band = conv(biofilm, kernel)[2:end-1, 2:end-1] .> 1e-15 # Some numbers are very small
    # return plotLogicalGrid(detachment_grid, Narrow_band) # Plot test for Narrow_band

    # Far is the inner part of the aggregate
    Far = biofilm .& .!Narrow_band
    # return plotLogicalGrid(detachment_grid, Far) # Plot test for Far

    # Inside of aggregate set to T infinity
    T[Far] .= Inf
    
    # Calculate detachment time for the narrow_band
    coordinates = findall(Narrow_band)
    for coor in coordinates
        y_index = coor[1]
        x_index = coor[2]
        Fdetach = calculateLocalDetachmentRate(x_index, y_index, kDet, grid, x_centre, y_centre)

        # COPIED IMPORTANT:
        # --------- IMPORTANT ----------
        # What is the impact of the number of free neighbours?
        # In idynomics it is calculated with the number of non-biomass
        # gridcells, but mathematically it does not make too much sense...?
        # --------END IMPORTANT --------

        nFreeNb = getFreeNeighbourCount(x_index, y_index, Visited)
        T[y_index, x_index] = grid.dx / (Fdetach * nFreeNb) # [h]
    end

    queue_size = sum(Narrow_band) # Number of items in queue
    T_queue = PriorityQueue()
    start_indices = findall(Narrow_band)
    # Add all items to the queue
    for index in start_indices
        push!(T_queue, string(index[1], ",", index[2]) => T[index])
    end

    # ----------------Fast Marching algorithm----------------------
    while queue_size > 0
        # Extract information of lowest T_value
        indices, value = first(T_queue)
        y_index, x_index = string_to_coordinates(indices)

        # Remove from Narrow_band and T_queue
        Narrow_band[y_index, x_index] = 0
        dequeue!(T_queue)
        queue_size = queue_size - 1

        # Add point to Visited
        Visited[y_index, x_index] = 1

        offSet = [0 0 1 -1; 1 -1 0 0]

        for nb in 1:length(offSet[1,:])
            x_index_nb = x_index + offSet[1, nb]
            y_index_nb = y_index + offSet[2, nb]
        
            if Visited[y_index_nb, x_index_nb] # If already visited, go to the next one
                continue
            end

            # recalculate T value
            T_val = recalculateT(T, x_index_nb, y_index_nb, kDet, grid, Visited, x_centre, y_centre)

            # If it was in Far, not in queue, so add to the queue
            if Far[y_index_nb, x_index_nb]
                Far[y_index_nb, x_index_nb] = 0 # remove from Far
                Narrow_band[y_index_nb, x_index_nb] = 1 # Add to Narrow_band
                push!(T_queue, string(y_index_nb, ",", x_index_nb) => T_val) # Add to the queue
                queue_size = queue_size + 1
            else    # if not in Far, it is in Narrow_band already, update value
                T_queue[string(y_index_nb, ",", x_index_nb)] = T_val
            end
            T[y_index_nb, x_index_nb] = T_val # Also add to the final storage
        end
    end

    # return plotLogicalGrid(detachment_grid, Visited) # Plot test for all Visited
    return plotLogicalGrid(detachment_grid, Far) # Plot test for Far


    return T
end

# ----------------------- START VISUALISATION -----------------------------
# # Read_file
# using Plots
# # Initilising 
# include(string(pwd(), "\\lib\\pre_processing\\create_mat.jl"))
# include(string(pwd(), "\\lib\\determine_where_bacteria_in_grid.jl"))
# filename = string(pwd(), "\\test\\test_file.xlsx")
# grid, bac, constants, settings, init_params = create_mat(filename)
# grid2bac, grid2nBacs = determine_where_bacteria_in_grid(grid, bac)
# calcTimeOfDetach(bac, grid, grid2bac, grid2nBacs, constants)

# ----------------------- END VISUALISATION -----------------------------