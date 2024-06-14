function plot_concentration_domain(sim_number, time_index, compound_indices)
    """
    """

    directory = string(pwd(), @sprintf("\\results\\%04d\\results2D.jld2", sim_number))

    bac_saved_vecint, bac_saved_matfloat , conc_saved  = load(directory, "bac_saved_vecint", "bac_saved_matfloat", "conc_saved")

    if time_index == -1
        time_index = findlast(bac_saved_vecint.nBacs .!= 0)
    end

    compound_strings = ["Butyrate", "Acetate", "CH4", "CO2", "H2", "NH3", "Na", "CH4_g", "CO2_g", "H2_g"]

    xmin = minimum(bac_saved_matfloat.x[time_index, 1:bac_saved_vecint.nBacs[time_index]]) * 1e6 / 2 - 30
    xmax = maximum(bac_saved_matfloat.x[time_index, 1:bac_saved_vecint.nBacs[time_index]]) * 1e6 / 2 + 30
    ymin = minimum(bac_saved_matfloat.y[time_index, 1:bac_saved_vecint.nBacs[time_index]]) * 1e6 / 2 - 30
    ymax = maximum(bac_saved_matfloat.y[time_index, 1:bac_saved_vecint.nBacs[time_index]]) * 1e6 / 2 + 30

    for i in eachindex(compound_indices)
        current_compound = compound_indices[i]
        heatmap(conc_saved[time_index,:,:,current_compound],
             label=compound_strings[current_compound],
             xlabel="Position along x-axis [µm]", 
             ylabel="Position along y-axis [µm]",
             xticks=(1:50:size(conc_saved[time_index,:,:,current_compound], 2), 1:100:(size(conc_saved[time_index,:,:,current_compound], 2) * 2)),
             yticks=(1:50:size(conc_saved[time_index,:,:,current_compound], 1), 1:100:(size(conc_saved[time_index,:,:,current_compound], 1) * 2)),
             xlim=(xmin,xmax),
             ylim=(ymin,ymax),
             title="$(compound_strings[current_compound]) at $((time_index-1) * 12) h",
             aspect_ratio=1)
        savefig(string(pwd(), "\\results\\$(sim_number)\\$(compound_strings[current_compound]) at $((time_index-1) * 12) h .png"))
    end

end
