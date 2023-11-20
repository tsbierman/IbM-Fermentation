function create_chunks(nChunks_dir, focus_region)
    """
    Creates chunks in the focus_region
    nChunks_dir is the number of Chunks per direction (amount of cores)
    focus_region is a rectangular region of interest in the entire domain

    chunks is a struct containing indices in x and y direction per chunk
    """

    dx = focus_region.x1 - focus_region.x0 + 1 # amount of grid cells
    dx_chunk = ceil(dx/nChunks_dir)            # grid cells per chunk x_direction

    dy = focus_region.y1 - focus_region.y0 + 1 # amount of grid cells
    dy_chunk = ceil(dy/nChunks_dir)            # grid cells per chunk y_direction

    # Calculate which start & end indices per chunk in x-direction
    indices_x = zeros(nChunks_dir, 2) # Store all starts and ends for x
    indices_y = zeros(nChunks_dir, 2) # Store all starts and ends for y

    # Store the beginning and ending of all the chunks
    let temp_x = focus_region.x0, temp_y = focus_region.y0
        for index_Chunks in eachindex(nChunks_dir)
            indices_x[index_Chunks, :] = [temp_x, temp_x + dx_chunk - 1]
            temp_x = temp_x + dx_chunk
            indices_y[index_Chunks, :] = [temp_y, temp_y + dy_chunk - 1]
            temp_y = temp_y + dy_chunk
        end
    end

    indices_x[end] = focus_region.x1
    indices_y[end] = focus_region.y1

    # recalibrate so that indices start at 1
    #   -> assumes that focus region is already being extracted and that
    #   chunks are only additional
    indices_x = indices_x .- focus_region.x0 .+ 1
    indices_y = indices_y .- focus_region.y0 .+ 1

    chunks = General()
    chunks.indices_x = indices_x
    chunks.indices_y = indices_y
    chunks.dx_chunk = dx_chunk
    chunks.dy_chunk = dy_chunk

    return chunks
end
