function sort_bacteria_into_chunks!(bac, grid, chunks, focus_region, nChunks_dir)
    """
    This function reorganises the indices in the bac for easier access in parallel computing
    Bacteria in the same chunk are next to each other

    Arguments
    bac:                A "General" struct containing all parameters related to the bacteria
    grid:               A "General" struct containing all parameters related to the grid
    chunks:             A "General" struct containing the start and end cooridinates of the chunks
    focus_region:       A "General" struct containing the x and y coordinates of the focus_region (diffusion region)
    nChunks_dir:        The amount of chunks that are generated in each direction
    
    Returns
    bac:                A reordened bac struct   
    """

    # Calculate which gridcell each bacteria is in
    ix = ceil.(bac.x / grid.dx)
    iy = ceil.(bac.y / grid.dy)

    # Calculate which chunk each bacteria is in
    ixChunk = ceil.((ix - focus_region.x0) / chunks.dx_chunk)
    iyChunk = ceil.((iy - focus_region.y0) / chunks.dy_chunk)

    # Combine chunk indices to a combined index
    bac_chunk = nChunks_dir * (ixChunk - 1) + iyChunk

    # sort chunks
    sortChunkIndex = sortperm(bac_chunk)

    # reorganise bac struct
    bac.x = bac.x[sortChunkIndex]
    bac.y = bac.y[sortChunkIndex]
    bac.species = bac.species[sortChunkIndex]
    bac.molarMass = bac.molarMass[sortChunkIndex]
    bac.radius = bac.radius[sortChunkIndex]
    bac.active = bac.active[sortChunkIndex]
    bac.mu = bac.mu[sortChunkIndex]

    return bac

end
