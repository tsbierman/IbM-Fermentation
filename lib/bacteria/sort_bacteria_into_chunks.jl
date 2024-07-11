function sort_bacteria_into_chunks!(bac_vecfloat, bac_vecint, bac_vecbool, grid_float, chunks_int, focus_region, nChunks_dir)
    """
    This function reorganises the bacteria (indices) in the structs for easier access in parallel computing.
    Bacteria in the same chunk are next to each other

    Arguments
    bac_XYZ:            A struct containing bacterial parameters
    grid_float:         A "Float" struct containing grid parameters of type Float
    chunks_int:         An "Int" struct containing the dx and dy of a chunk
    focus_region:       An "Int" struct containing the x and y coordinates of the focus_region (diffusion region)
    nChunks_dir:        An Integer, indicating the amount of chunks that are generated in each direction
    
    Returns
    bac_XYZ:            A struct containing bacterial parameters
    """

    # Calculate in which gridcell each bacteria is
    ix = ceil.(bac_vecfloat.x ./ grid_float.dx)
    iy = ceil.(bac_vecfloat.y ./ grid_float.dy)

    # Calculate in which chunk each bacteria is
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
    bac_vecint.colony_nums = bac_vecint.colony_nums[sortChunkIndex]

    return bac_vecfloat, bac_vecint, bac_vecbool
end
