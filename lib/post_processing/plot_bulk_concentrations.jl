function plot_bulk_concentrations(sim_number, compound_indices)
    """ 
    This function plots bulk concentration for a single simulation. Each compound gets its own graph

    Arguments
    sim_number:             The simulation number
    coumpound_indecies:     A vector or the indices that need to be plotted (based on how they are stored)

    Returns
    A lineplot (and saves it)
    """

    # Get data
    if isfile(string(pwd(), @sprintf("\\results\\%04d\\results1D.jld2", sim_number)))
        directory = string(pwd(), @sprintf("\\results\\%04d\\results1D.jld2", sim_number))
    else
        directory = string(pwd(), @sprintf("\\results\\%04d\\results2D.jld2", sim_number))
    end

    bac_saved_vecint, reactor_saved_matfloat  = load(directory, "bac_saved_vecint", "reactor_saved_matfloat")

    lastnonzero = findlast(bac_saved_vecint.nBacs .!= 0)
    bulk_conc_to_plot = reactor_saved_matfloat.bulk_concs[1:lastnonzero, :]

    # Data has been stored every 12 hours
    plot_time = (0:lastnonzero-1) * 12

    # Compound names for Anaerobic Fermentation
    compound_strings = ["Butyrate", "Acetate", "CH4", "CO2", "H2", "NH3", "Na", "CH4_g", "CO2_g", "H2_g"]

    for i in eachindex(compound_indices)
        current_compound = compound_indices[i]
        plot(plot_time, bulk_conc_to_plot[1:lastnonzero, current_compound], label=compound_strings[current_compound],
             xlabel="Time (h)",
             ylabel="Concentration (M)",
             title="$(compound_strings[current_compound])",
             xticks=(0:24:maximum(plot_time), 0:24:maximum(plot_time)), # Ticks on x-axis every 24 hours
             ylim=(0, maximum(bulk_conc_to_plot[1:lastnonzero, current_compound])*1.05))
        savefig(string(pwd(), "\\results\\$(sim_number)\\$(compound_strings[current_compound]).png"))
    end
end

function plot_bulk_concentrations_multiple(sim_numbers, compound_indices)
    """ 
    This function plots bulk concentrations for multiple simulation in the same figure.
    All compounds end up in the same figure.

    Arguments
    sim_numbers:                A vector with the simulation number_of_bacs
    compound_indices:           A vector with the compound indices (based on how they are stored)

    Returns
    A lineplot (and saves it)
    """

    # Initialise clean plot
    plot()

    # Compound Names for Anaerobic Fermentation
    compound_strings = ["Butyrate", "Acetate", "CH4", "CO2", "H2", "NH3", "Na", "CH4_g", "CO2_g", "H2_g"]

    # Colors and linestyles used in Anaerobic Fermentation
    colors = ["#276419", "#DE77AE", "#B8E186"]
    linestyles = [:solid, :dash]
    labels = ["Rosenbrock23", "Rodas4P"]
    
    # Plot legend a single time
    for i in eachindex(compound_indices)
        plot!([],[], color = colors[2], label=compound_strings[compound_indices[i]])
    end

    for j in eachindex(sim_numbers)
        sim_number = sim_numbers[j]

        # Get data
        if isfile(string(pwd(), @sprintf("\\results\\%04d\\results1D.jld2", sim_number)))
            directory = string(pwd(), @sprintf("\\results\\%04d\\results1D.jld2", sim_number))
        else
            directory = string(pwd(), @sprintf("\\results\\%04d\\results2D.jld2", sim_number))
        end

        bac_saved_vecint, reactor_saved_matfloat  = load(directory, "bac_saved_vecint", "reactor_saved_matfloat")

        lastnonzero = findlast(bac_saved_vecint.nBacs .!= 0)
        bulk_conc_to_plot = reactor_saved_matfloat.bulk_concs[1:lastnonzero, :] .* 1000 # Convert to mM

        plot_time = (0:lastnonzero-1) * 12

        for i in eachindex(compound_indices)
            current_compound = compound_indices[i]
            plot!(plot_time, bulk_conc_to_plot[1:lastnonzero, current_compound],
                #  label=compound_strings[current_compound],
                xlabel="Time (h)",
                ylabel="Concentration (mM)",
                linestyle = linestyles[j],
                label=nothing,
                legend_position =:best,
                #  title="$(compound_strings[current_compound])",
                c = colors[2],
                # xlim=(0, 50),
                xticks=(0:240:2400, 0:240:2400), # Ticks every 240 hours
                ylim=(0, maximum(bulk_conc_to_plot[1:lastnonzero, current_compound])*1.05))

            # savefig(string(pwd(), "\\results\\$(sim_number)\\$(compound_strings[current_compound]).png"))
        end
    end

    # The following code can be used to run the ADM1 and plot it together with the IbM results in a single graph
    # ADM1_sol = run_ADM1()
    # time, but_conc, ace_conc, CH4_conc, CO2_conc, CH4_g_conc, CO2_g_conc = ADM1_sol[1], ADM1_sol[2], ADM1_sol[3], ADM1_sol[5], ADM1_sol[6], ADM1_sol[13], ADM1_sol[14]
    # plot!(time, but_conc .* 1000, linestyle=:dash, c=colors[1], label=nothing)
    # plot!(time, ace_conc .* 1000, linestyle=:dash, c=colors[3], label=nothing)
    # plot!(time, CO2_conc .* 1000, linestyle=:dash, c=colors[2], label=nothing)
    # plot!(time, CH4_conc .* 1000, linestyle=:dot, c=colors[2], label=nothing)
    # plot!(time, CH4_g_conc .* 1000, linestyle=:dot, c=colors[2], label=nothing)
    # plot!(time, CO2_g_conc .* 1000, linestyle=:dot, c=colors[2], label=nothing)

    savefig(string(pwd(), "\\results\\$(sim_numbers).png"))
end
