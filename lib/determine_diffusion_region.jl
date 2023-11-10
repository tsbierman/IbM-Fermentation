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

    # Continue here

    return diffusionNodesX, diffusionNodesY
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

    return diffusion_region, focus_region
end