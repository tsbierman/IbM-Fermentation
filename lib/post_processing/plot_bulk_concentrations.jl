function plot_bulk_concentrations(sim_number, compound_indices)
    """ 
    """

    if isfile(string(pwd(), @sprintf("\\results\\%04d\\results1D.jld2", sim_number)))
        directory = string(pwd(), @sprintf("\\results\\%04d\\results1D.jld2", sim_number))
    else
        directory = string(pwd(), @sprintf("\\results\\%04d\\results2D.jld2", sim_number))
    end

    bac_saved_vecint, reactor_saved_matfloat  = load(directory, "bac_saved_vecint", "reactor_saved_matfloat")

    lastnonzero = findlast(bac_saved_vecint.nBacs .!= 0)
    bulk_conc_to_plot = reactor_saved_matfloat.bulk_concs[1:lastnonzero, :]

    plot_time = (0:lastnonzero-1) * 12

    compound_strings = ["Butyrate", "Acetate", "CH4", "CO2", "H2", "NH3", "Na", "CH4_g", "CO2_g", "H2_g"]

    for i in eachindex(compound_indices)
        current_compound = compound_indices[i]
        plot(plot_time, bulk_conc_to_plot[1:lastnonzero, current_compound], label=compound_strings[current_compound],
             xlabel="Time (h)",
             ylabel="Concentration (M)",
             title="$(compound_strings[current_compound])",
             xticks=(0:24:maximum(plot_time), 0:24:maximum(plot_time)),
             ylim=(0, maximum(bulk_conc_to_plot[1:lastnonzero, current_compound])*1.05))
        savefig(string(pwd(), "\\results\\$(sim_number)\\$(compound_strings[current_compound]).png"))
    end
end
