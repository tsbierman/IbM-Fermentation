function getDiffusionNodes(bac, grid_float, grid_int)
    """
    This function determines in both X and Y direction which nodes are potentially in the diffusion region

    Arguments
    bac:                A "General" struct containing all parameters related to the bacteria
    grid:               A "General" struct containing all parameters related to the grid

    Returns
    diffusionNodesX:    A (nx,) vector indicating which columns are potentially in the diffusion region
    diffusionNodesY:    A (ny,) vector indicating which rows are potentially in the diffusion region
    """

    # Extract variables
    x_min = minimum(bac.x)
    x_max = maximum(bac.x)
    y_min = minimum(bac.y)
    y_max = maximum(bac.y)

    # Define an offset, so we guarantee that it will contain the correct indices
    offsetX = grid_float.blayer_thickness + grid_float.dx
    offsetY = grid_float.blayer_thickness + grid_float.dy

    nodeEndCoordinatesX = (1:grid_int.nx) * grid_float.dx
    nodeEndCoordinatesY = (1:grid_int.ny) * grid_float.dy

    # The following part is designed such that even if a small part of a grid cell
    # is within the diffusion region (boundary_layer + some offset), it is true
    diffusionNodesX = (nodeEndCoordinatesX .> x_min - offsetX) .& (nodeEndCoordinatesX .- grid_float.dx .< x_max + offsetX)
    diffusionNodesY = (nodeEndCoordinatesY .> y_min - offsetY) .& (nodeEndCoordinatesY .- grid_float.dy .< y_max + offsetY)

    return diffusionNodesX, diffusionNodesY
end

function isWithinBoundaryLayer(bac_x, bac_y, gridcell_centre, blayer_thickness)
    """
    This function determines whether a gridcell is within the boundary layer of a bacterium. 
    For this, only the centre of the gridcell is considered

    Arguments
    bac_x, bac_y:           The coordinates of the bacterium
    gridcell_centre:        The coordinates of the gridcell centre (x,y)
    blayer_thickness:       The boundary layer blayer_thickness [m]

    Returns
    isBlayer:               A boolean indicating whether the gridcell is within the boundary layer
    """

    isBlayer = blayer_thickness^2 > sum(([bac_x, bac_y] - gridcell_centre).^2)
    return isBlayer
end


function determine_focus_region(diffRegion)
    """
    This function determines the extraction (focus) region, which:
        - includes all diffusion nodes
        - has an odd number of nodes in x & y direction (for efficient diffusion)
        - has at least one bulk layer node at each side
    
    Arguments
    diffRegion:                 A BitMatrix indicating per gridcell whether that cell is in the diffusion region

    Returns
    extraction_region:          A "General" struct containing the start and end coordinates of the region of interest
    """

    # Struct to save
    extraction_region = Int_struct()

    # Determine first and last x
    x_diffRegion = sum(diffRegion, dims= 1)
    first_x = findfirst(x_diffRegion .!= 0)[2] - 1
    last_x = findlast(x_diffRegion .!= 0)[2] + 1

    dx = last_x - first_x + 1           # number of diffusion cells
    if mod(dx, 2) == 0
        first_x = first_x - 1           # make uneven
    end

    # Determine first and last y
    y_diffRegion = sum(diffRegion, dims =2)
    first_y = findfirst(y_diffRegion .!= 0)[1] - 1
    last_y = findlast(y_diffRegion .!= 0)[1] + 1

    dy = last_y - first_y + 1           # number of diffusion cells
    if mod(dy, 2) == 0
        first_y = first_y - 1           # make uneven
    end

    if first_x < 1 || last_x > size(diffRegion, 2) || first_y < 1 || last_y > size(diffRegion, 1)
        @warn "DEBUG:actionRequired, debug: not enough bulk liquid present around granule"
    end

    extraction_region.x0 = first_x
    extraction_region.x1 = last_x
    extraction_region.y0 = first_y
    extraction_region.y1 = last_y

    # DEBUG
    # if any(diffRegion[:, [1,end]] .== 1) || any(diffRegion[[1,end], :] .== 1)
    #     xRange = first_x:last_x
    #     yRange = first_y:last_y
    #     display(diffRegion[yRange, xRange])
    # end
    # END DEBUG

    return extraction_region
end


function determine_diffusion_region(grid2bac, grid2nBacs, bac, grid_float, grid_int)
    """
    This function determines which gridcells are in the diffusion layer. The diffusion layer is defined
    as the gridcells that are within boundary layer distance from a bacterium.

    Arguments
    grid2bac:           A matrix (ny, nx, ?) which contains for each gridcell which bacteria is located
                        there. The number corresponds to the index in the bac struct
    grid2nBacs:         A (ny, nx) matrix which contains for each gridcell how many bacteria are located there
    bac:                A "General" struct containing all parameters related to the bacteria
    grid:               A "General" struct containing all parameters related to the grid
    
    Returns
    diffusion_region:   A BitMatrix indicating per gridcell whether that cell is in the diffusion region
    focus_region:       A "General" struct containing the start and end coordinates of the region of interest
    """
    
    # Initialise
    diffusion_region = BitArray(zeros(grid_int.ny, grid_int.nx))

    # Get diffusion region
    diffNodesX, diffNodesY = getDiffusionNodes(bac, grid_float, grid_int)

    # In the diffusion region, apply convolution to find the boundary of gridcells with bacteria
    kernel = ones(3,3) ./ -8
    kernel[2,2] = 1
    hasBac = grid2nBacs[diffNodesY, diffNodesX] .> 0                # A BitArray of only estimated diffusion region
    isBacBoundary = conv(hasBac, kernel)[2:end-1,2:end-1] .> 1e-15  # Some numbers end up very small, but technically larger than 0 (1e-17)
    
    # For the boundary of bacterial grid cells, compute for the neighbouring grid cells whether
    # they are in diffusion region
    maxOffsetX = ceil(grid_float.blayer_thickness / grid_float.dx)              # Maximum offset to check
    dx = findfirst(diffNodesX) - 1                                              # Adjust for difference of total and diffusion grid
    maxOffsetY = ceil(grid_float.blayer_thickness / grid_float.dy)              # Maximum offset to check
    dy = findfirst(diffNodesY) - 1                                              # Adjust for difference of total and diffusion grid

    cart_indices = findall(isBacBoundary)                                       # Array of Cartesian indices
    isDiffRegion = deepcopy(hasBac)

    for index in cart_indices                                       # For every boundary cell
        for di in -maxOffsetX:maxOffsetX                            # Check X-direction
            for dj in -maxOffsetY:maxOffsetY                        # Check Y-direction
                if isDiffRegion[Int(index[1] + dj), Int(index[2] + di)] != 0  # Check whether already in diffusion
                    continue                                        # Skip if already in diffusion
                else
                    # perform actual check
                    bac_indices = findall(grid2bac[index[1] + dy, index[2] + dx, :] .!= 0)  # Find how many bacteria spots are occupied in the grid2bac and by which
                    bacs = grid2bac[index[1] + dy, index[2] + dx, bac_indices]              # Get bacterial indices (greater picture)
                    gridcell_centre = [grid_float.dx * (dx + di + index[2]) - grid_float.dx/2, grid_float.dy * (dy + dj + index[1]) - grid_float.dy/2]  # Calculate centre of desired grid

                    for iBac in bacs                                # For every of the bacteria
                        if isWithinBoundaryLayer(bac.x[iBac], bac.y[iBac], gridcell_centre, grid_float.blayer_thickness)              # Check whether the grid cell is within reach
                            isDiffRegion[Int(index[1] + dj), Int(index[2] + di)] = 1
                            break
                        end
                    end
                end
            end
        end
    end

    # Place back in original matrix
    diffusion_region[diffNodesY, diffNodesX] = isDiffRegion

    # Determine focus_region
    focus_region = determine_focus_region(diffusion_region)

    return diffusion_region, focus_region
end
