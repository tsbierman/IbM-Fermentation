function getDiffusionNodes(bac, grid)
    """
    Determine in both X and Y direction which nodes are potentially in the diffusion region

    bac is a struct containing all information regarding the bacteria
    grid is a struct containing all information regarding the grid
    It returns 2 logical arrays that signify per direction which nodes are potentially in the diffusion region
    """

    # Extract variables
    x_min = minimum(bac.x)
    x_max = maximum(bac.x)
    y_min = minimum(bac.y)
    y_max = maximum(bac.y)

    # Define an offset, so we guarantee that it will contain the 
    offsetX = grid.blayer_thickness + grid.dx
    offsetY = grid.blayer_thickness + grid.dy

    nodeEndCoordinatesX = (1:grid.nx) * grid.dx
    nodeEndCoordinatesY = (1:grid.ny) * grid.dy

    # The following part is designed such that even if a small part of a grid cell
    # is within the diffusion region (boundary_layer + some offset), it is taken into account
    diffusionNodesX = (nodeEndCoordinatesX .> x_min - offsetX) .& (nodeEndCoordinatesX .- grid.dx .< x_max + offsetX)
    diffusionNodesY = (nodeEndCoordinatesY .> y_min - offsetY) .& (nodeEndCoordinatesY .- grid.dy .< y_max + offsetY)

    return diffusionNodesX, diffusionNodesY
end

function isWithinBoundaryLayer(bac_x, bac_y, gridcell_centre, blayer_thickness)
    """
    Determine whether cell is within boundary layer. For this, only the centre is considered
    bac_x, bac_y are bacterial coordinates
    gridcell_centre is centre coordinates (x,y)
    blayer_thickness is the boundary layer blayer_thickness
    returns a boolean
    """

    # isBlayer = sqrt((bac_x - gridcell_centre[1])^2 + (bac_y - gridcell_centre[2])^2) <= blayer_thickness
    isBlayer = blayer_thickness^2 > sum(([bac_x, bac_y] - gridcell_centre).^2)
    return isBlayer
end


function determine_focus_region(diffRegion)
    """
    Determine extraction region, which:
    - includes all diffusion nodes
    - has an odd number of nodes in x & y direction (efficient diffusion)
    - has at least one bulk layer node at each side
    """
    # Struct to save
    extraction_region = General()

    # Determine first and last x
    x_diffRegion = sum(diffRegion, dims= 1)
    first_x = findfirst(x_diffRegion .!= 0)[2] - 1
    last_x = findlast(x_diffRegion .!= 0)[2] + 1

    dx = last_x - first_x + 1 # number of diffusion cells
    if mod(dx, 2) == 0
        first_x = first_x - 1 # make uneven
    end

    # Determine first and last y
    y_diffRegion = sum(diffRegion, dims =2)
    first_y = findfirst(y_diffRegion .!= 0)[1] - 1
    last_y = findlast(y_diffRegion .!= 0)[1] + 1

    dy = last_y - first_y + 1 # number of diffusion cells
    if mod(dy, 2) == 0
        first_y = first_y - 1 # make uneven
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


function determine_diffusion_region(grid2bac, grid2nBacs, bac, grid)
    """
    Determine which grid cells are in the diffusion layer. The diffusion layer is defined
    as the grid cells that are within boundary distance from a bacterium.

    grid2bac is a nx * ny * ? matrix with the index of every bacterium in each grid cell
    grid2nBacs is a nx * ny matrix with the number of bacteria in the grid cell
    bac is a struct containing all information regarding the bacteria
    grid is a struct containing all information regarding the grid
    
    Diffusion region is a logical matrix which indicates which grid cells are within diffusion region
    """
    
    # Initialise
    diffusion_region = BitArray(zeros(grid.ny, grid.nx))

    # Get diffusion region
    diffNodesX, diffNodesY = getDiffusionNodes(bac, grid)

    # In the diffusion region, appy convolution to find boundary of grid with cells with bacteria
    kernel = ones(3,3) ./ -8
    kernel[2,2] = 1
    hasBac = grid2nBacs[diffNodesY, diffNodesX] .> 0 # BitArray of only estimated diffusion region
    isBacBoundary = conv(hasBac, kernel)[2:end-1,2:end-1] .> 1e-15 #some numbers end up very small, but technically larger than 0 (1e-17)
    
    # For the boundary of bacterial grid cells, compute for the neighbouring grid cells whether
    # they are in diffusion region
    maxOffsetX = ceil(grid.blayer_thickness / grid.dx)      # Maximum offset to check
    dx = findfirst(diffNodesX) - 1                          # Adjust for difference of total and diffusion grid
    maxOffsetY = ceil(grid.blayer_thickness / grid.dy)      # Maximum offset to check
    dy = findfirst(diffNodesY) - 1                          # Adjust for difference of total and diffusion grid

    cart_indices = findall(isBacBoundary) # Array of Cartesian indices
    isDiffRegion = deepcopy(hasBac)

    for index in cart_indices                                       # For every boundary cell
        for di in -maxOffsetX:maxOffsetX                            # Check X-direction
            for dj in -maxOffsetY:maxOffsetY                        # Check Y-direction
                if isDiffRegion[Int(index[1] + dj), Int(index[2] + di)] != 0  # Check whether already in diffusion
                    continue                                        # Skip if already in diffusion
                else
                    # perform actual check
                    bac_indices = findall(grid2bac[index[1] + dy, index[2] + dx, :] .!= 0) # Find which bacteria spots are occupied in the grid cell
                    bacs = grid2bac[index[1] + dy, index[2] + dx, bac_indices]           # Get bacterial indices (greater picture)
                    gridcell_centre = [grid.dx * (dx + di + index[2]) - grid.dx/2, grid.dy * (dy + dj + index[1]) - grid.dy/2] # Calculate centre of desired grid

                    for iBac in bacs    # For every of the bacteria
                        if isWithinBoundaryLayer(bac.x[iBac], bac.y[iBac], gridcell_centre, grid.blayer_thickness) # Check whether the grid cell is within reach
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
