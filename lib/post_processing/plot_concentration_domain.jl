function plot_concentration_domain(sim_number, time_index, compound_indices)
    """
    This function plots a heatmap of the specified compounds. Each compound is plotted in their own figure.
    This function only works when a 2-Dimensional map is stored during the simulation.

    Arguments:
    sim_number:         The simulation number
    time_index:         The index of the stored data to be plotted. If time_index = -1, it uses the most recent timepoint.
    compound_indices:   A vector of the indices of the compounds that need to be plotted (based on the order they are saved in)

    Returns:
    A heatmap (and saves it)
    """

    # Get data
    directory = string(pwd(), @sprintf("\\results\\%04d\\results2D.jld2", sim_number))

    bac_saved_vecint, bac_saved_matfloat , conc_saved  = load(directory, "bac_saved_vecint", "bac_saved_matfloat", "conc_saved")

    # Get most recent data
    if time_index == -1
        time_index = findlast(bac_saved_vecint.nBacs .!= 0)
    end

    # Compound names for Anaerobic Fermentation
    compound_strings = ["Butyrate", "Acetate", "CH4", "CO2", "H2", "NH3", "Na", "CH4_g", "CO2_g", "H2_g"]

    # Set boundaries for plotting (in µm)
    xmin = minimum(bac_saved_matfloat.x[time_index, 1:bac_saved_vecint.nBacs[time_index]]) * 1e6 / 2 - 30
    xmax = maximum(bac_saved_matfloat.x[time_index, 1:bac_saved_vecint.nBacs[time_index]]) * 1e6 / 2 + 30
    ymin = minimum(bac_saved_matfloat.y[time_index, 1:bac_saved_vecint.nBacs[time_index]]) * 1e6 / 2 - 30
    ymax = maximum(bac_saved_matfloat.y[time_index, 1:bac_saved_vecint.nBacs[time_index]]) * 1e6 / 2 + 30

    for i in eachindex(compound_indices)
        current_compound = compound_indices[i]
        heatmap(conc_saved[time_index,:,:,current_compound] .* 1000, # Plotting in mM
             label=compound_strings[current_compound],
             xlabel="Position along x-axis [µm]", 
             ylabel="Position along y-axis [µm]",
             colorbar_title= " \n Butyrate concentration (mM)", # Extra newline to prevent overlap
             right_margin = 2Plots.mm,
             xticks=(1:50:size(conc_saved[time_index,:,:,current_compound], 2), 1:100:(size(conc_saved[time_index,:,:,current_compound], 2) * 2)),
             yticks=(1:50:size(conc_saved[time_index,:,:,current_compound], 1), 1:100:(size(conc_saved[time_index,:,:,current_compound], 1) * 2)),
             xlim=(xmin,xmax),
             ylim=(ymin,ymax),
            #  clim=(0,1.13),
            #  title="$(compound_strings[current_compound]) at $((time_index-1) * 12) h",
             aspect_ratio=1)
        savefig(string(pwd(), "\\results\\$(sim_number)\\$(compound_strings[current_compound]) at $((time_index-1) * 12) h .png"))
    end
end
