function sort_bacteria_into_chunks!(bac_vecfloat, bac_vecint, bac_vecbool, grid_float, chunks_int, focus_region, nChunks_dir)
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
    ix = ceil.(bac_vecfloat.x ./ grid_float.dx)
    iy = ceil.(bac_vecfloat.y ./ grid_float.dy)

    # Calculate which chunk each bacteria is in
    ixChunk = ceil.((ix - focus_region.x0) / chunks_int.dx_chunk)
    iyChunk = ceil.((iy - focus_region.y0) / chunks_int.dy_chunk)

    # Combine chunk indices to a combined index
    bac_chunk = nChunks_dir * (ixChunk - 1) + iyChunk

    # sort chunks
    sortChunkIndex = sortperm(bac_chunk)

    # reorganise bac struct
    bac_vecfloat.x = bac_vecfloat.x[sortChunkIndex]
    bac_vecfloat.y = bac_vecfloat.y[sortChunkIndex]
    bac_vecint.species = bac_vecint.species[sortChunkIndex]
    bac_vecfloat.molarMass = bac_vecfloat.molarMass[sortChunkIndex]
    bac_vecfloat.radius = bac_vecfloat.radius[sortChunkIndex]
    bac_vecbool.active = bac_vecbool.active[sortChunkIndex]
    bac_vecfloat.mu = bac_vecfloat.mu[sortChunkIndex]

    return bac_vecfloat, bac_vecint, bac_vecbool

end
