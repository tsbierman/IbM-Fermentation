function plot_growth_rate_heatmap(sim_number, time_index, finished, specie)
    """
    This function generates a granule with on top of it a heatmap of the growth rate of a specified specie.

    Arguments
    sim_number:             The simulation number
    time_index:             The index of which data should be used. If time_index = -1, it uses the most recent data
    finished:               A boolean indication whether the simulation has been completed
    specie:                 An Integer indicating of which specie the growth rate should be plotted

    Returns
    An aggregate with growth rate heat map
    """

    # Get directory and simulation_file
    directory = string(pwd(), @sprintf("\\results\\%04d\\results2D.jld2", sim_number))

    if finished
        sim_file = string(pwd(), @sprintf("\\results\\%04d\\sim_%04d.jld2", sim_number, sim_number))
    else
        sim_file = string(pwd(), @sprintf("\\sim_%04d.jld2", sim_number))
    end

    # Extract grid, bacteria and concentration information
    grid_float, grid_int = load(sim_file, "grid_float", "grid_int")
    bac_saved_vecint, bac_saved_matint, bac_saved_matfloat = load(directory, "bac_saved_vecint", "bac_saved_matint", "bac_saved_matfloat")

    if time_index == -1
        time_index = findlast(bac_saved_vecint.nBacs .!= 0)
    end

    # Find most recent index
    last_non_zero = findlast(bac_saved_vecint.nBacs .!= 0)

    # Select the data
    nBacs = bac_saved_vecint.nBacs[time_index]
    species = bac_saved_matint.species[time_index, 1:nBacs]

    short_bac_vecfloat = VectorFloat_struct()
    bac_x = bac_saved_matfloat.x[time_index, 1:nBacs]
    bac_y = bac_saved_matfloat.y[time_index, 1:nBacs]
    mu_s = bac_saved_matfloat.mu[time_index, 1:nBacs]
    short_bac_vecfloat.x = bac_x
    short_bac_vecfloat.y = bac_y

    # Determine bacterial position
    grid2bac, grid2nBacs = determine_where_bacteria_in_grid(grid_float, grid_int, short_bac_vecfloat)

    # Initialise storage
    growth_array = zeros(size(grid2nBacs, 1) * 2, size(grid2nBacs, 2) * 2)
    growth_array[:,:] .= NaN

    # Traverse over all grid cells
    for ix in axes(grid2nBacs, 2)
        for iy in axes(grid2nBacs, 1)

            if grid2nBacs[iy, ix] > 0 # Only when bacteria are present
                number_of_bacs = grid2nBacs[iy, ix]

                for bac_num in 1:number_of_bacs
                    # Get bac index in whole system, its specie number and the catabolism
                    bac_index = grid2bac[iy,ix,bac_num]
                    bac_specie = species[bac_index]

                    if bac_specie == specie # If correct specie, store growth rate
                        growth_array[iy*2-1:iy*2,ix*2-1:ix*2] .= mu_s[bac_index]
                    end
                end
            end
        end
    end

    # Determine plotting range
    xmin = minimum(bac_saved_matfloat.x[last_non_zero, 1:bac_saved_vecint.nBacs[time_index]]) * 1e6 - 10
    xmax = maximum(bac_saved_matfloat.x[last_non_zero, 1:bac_saved_vecint.nBacs[time_index]]) * 1e6 + 10
    ymin = minimum(bac_saved_matfloat.y[last_non_zero, 1:bac_saved_vecint.nBacs[time_index]]) * 1e6 - 10
    ymax = maximum(bac_saved_matfloat.y[last_non_zero, 1:bac_saved_vecint.nBacs[time_index]]) * 1e6 + 10

    # Plot aggregate
    plotBacs(sim_number, finished, time_index, false)

    # Add heatmap
    heatmap!(growth_array,
            xlim=(xmin,xmax),
            ylim=(ymin,ymax),
            right_margin = 2Plots.mm,
            colorbar_title= " \nGrowth rate (h-1)",
            # clim=(-0.0008, 0.0023),
            aspect_ratio=1)

    savefig(string(pwd(), "\\results\\mu_heat_$(sim_number)_$(time_index)_$(specie).png"))
end
