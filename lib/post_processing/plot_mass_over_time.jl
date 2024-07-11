function plot_mass_over_time(sim_number)
    """
    This function plots the active mass (in moles) of the organisms over time. Additionally, it plots their fractions.

    Arguments
    sim_number:             The simulation number

    Returns
    Two lineplots (and saves them)
    """

    # Get data
    if isfile(string(pwd(), @sprintf("\\results\\%04d\\results1D.jld2", sim_number)))
        directory = string(pwd(), @sprintf("\\results\\%04d\\results1D.jld2", sim_number))
    else
        directory = string(pwd(), @sprintf("\\results\\%04d\\results2D.jld2", sim_number))
    end

    bac_vecint, bac_matfloat, bac_matint, bac_matbool = load(directory, "bac_saved_vecint", "bac_saved_matfloat", "bac_saved_matint", "bac_saved_matbool")

    # Get last time point
    lastnonzero = findlast(bac_vecint.nBacs .!= 0)

    # Bacterial names and colors for Anaerobic Fermentation
    bac_names = ["BO" "AM" "HM"]
    color_list = ["#A020F0" "#E8BB00"  "#3bd6c6"]
    
    # Create storage
    store_fractions = zeros(size(bac_names, 2), lastnonzero)
    store_mass = zeros(size(bac_names, 2), lastnonzero)

    for j in 1:lastnonzero
        # Select data
        nbac = bac_vecint.nBacs[j]
        actives = bac_matbool.active[j,1:nbac]
        species = bac_matint.species[j,1:nbac]
        radius = bac_matfloat.radius[j,1:nbac]

        # Calculate moles and fractions and store them
        total_moles_per_specie = [sum((radius[actives .& (species .== i)]).^3 .* (4 ./ 3) .* pi) for i in 1:maximum(species)] .* 500000 ./ 22.6232 
        total_moles = sum(total_moles_per_specie)
        mole_fractions = total_moles_per_specie ./ total_moles
        store_mass[:,j] = total_moles_per_specie
        store_fractions[:,j] = mole_fractions
    end

    plot_time = (0:lastnonzero-1) * 12

    # Create clean plot
    plot()

    # Plot the mass of each specie
    for bac in eachindex(bac_names)
        plot!(plot_time, store_mass[bac,:], color = color_list[bac], label=bac_names[bac],
        xlabel="Time (h)", ylabel="Mass (mole)", xticks=(0:24:maximum(plot_time), 0:24:maximum(plot_time)))
    end
    savefig(string(pwd(), "\\results\\$(sim_number)\\$(sim_number)_mass_over_time.png"))

    # Create clean plot
    plot()

    # Plot the fraction of each specie
    for bac in eachindex(bac_names)
        plot!(plot_time, store_fractions[bac,:], color = color_list[bac], label=bac_names[bac], 
            xlabel="Time (h)", ylabel="Fraction", xticks=(0:24:maximum(plot_time), 0:24:maximum(plot_time)), ylim=(0,1), xlim=(0,maximum(plot_time)),
            yticks=(0:0.1:1, 0:0.1:1))
    end
    savefig(string(pwd(), "\\results\\$(sim_number)\\$(sim_number)_fractions_over_time.png"))
end

function plot_mass_over_time_multiple(sim_numbers)
    """
    This function plots the active mass (in moles) of the organisms over time. Additionally, it plots their fractions.
    It plots the several simulation in a single graph with different linestyles.

    Arguments
    sim_numbers:             A vector with the simulation numbers

    Returns
    Two lineplots (and saves them)
    """

    # Names, colors and linestyles used in Anaerobic Fermentation
    bac_names = ["BO" "AM" "HM"]
    color_list = ["#A020F0" "#E8BB00" "#3bd6c6"]
    styles = [:solid, :dash]

    # Create new plot and generate legend only once
    plot()
    for i in eachindex(bac_names)
        plot!([],[], color = color_list[i], label=bac_names[i])
    end

    for sims in eachindex(sim_numbers)
        # Get data from saved files
        sim_number = sim_numbers[sims]
        if isfile(string(pwd(), @sprintf("\\results\\%04d\\results1D.jld2", sim_number)))
            directory = string(pwd(), @sprintf("\\results\\%04d\\results1D.jld2", sim_number))
        else
            directory = string(pwd(), @sprintf("\\results\\%04d\\results2D.jld2", sim_number))
        end

        bac_vecint, bac_matfloat, bac_matint, bac_matbool = load(directory, "bac_saved_vecint", "bac_saved_matfloat", "bac_saved_matint", "bac_saved_matbool")

        lastnonzero = findlast(bac_vecint.nBacs .!= 0)
        
        # Create storage
        store_fractions = zeros(size(bac_names, 2), lastnonzero)
        store_mass = zeros(size(bac_names, 2), lastnonzero)

        for j in 1:lastnonzero
            # Select data
            nbac = bac_vecint.nBacs[j]
            actives = bac_matbool.active[j,1:nbac]
            species = bac_matint.species[j,1:nbac]
            radius = bac_matfloat.radius[j,1:nbac]

            # Calculate moles and fractions
            total_moles_per_specie = [sum((radius[actives .& (species .== i)]).^3 .* (4 ./ 3) .* pi) for i in 1:maximum(species)] .* 500000 ./ 22.6232 
            total_moles = sum(total_moles_per_specie)
            mole_fractions = total_moles_per_specie ./ total_moles
            store_mass[:,j] = total_moles_per_specie                        # For moles
            # store_mass[:,j] = total_moles_per_specie ./ (5.75e-7) .* 1000 # For concentrations in mM
            store_fractions[:,j] = mole_fractions
        end

        plot_time = (0:lastnonzero-1) * 12

        # Plot the active masses
        for bac in eachindex(bac_names)
            plot!(plot_time, store_mass[bac,:], color = color_list[bac], label=nothing, xlim=(0,maximum(plot_time)), ylim=(0,maximum(store_mass)),
            xlabel="Time (h)", ylabel="Mass (mole)", xticks=(0:240:maximum(plot_time), 0:240:maximum(plot_time)), linestyle=styles[sims])
        end
        savefig(string(pwd(), "\\results\\$(sim_numbers)_mass_over_time.png"))


    # The following code can be used to plot the biomass concentrations of the ADM1 in the same plot. In order to use this, disable the fraction plotting (lines 148-153)
    # and move the savefig (line 134) to after the ADM1 plotting (after line 145).

    # sol = run_ADM1()
    # time, Xbut, Xac, Xh2 = sol[1], sol[8], sol[9], sol[10]

    # plot!(time, Xbut .* 1000, color=color_list[1], label=nothing, linestyle=:dash)
    # plot!(time, Xac .* 1000, color=color_list[2], label=nothing, linestyle=:dash)
    # plot!(time, Xh2 .* 1000, color=color_list[3], label=nothing, linestyle=:dash)

        # Plot the fractions
        for bac in eachindex(bac_names)
            plot!(plot_time, store_fractions[bac,:], color = color_list[bac], label=bac_names[bac], 
                xlabel="Time (h)", ylabel="Fraction", xticks=(0:24:maximum(plot_time), 0:24:maximum(plot_time)), ylim=(0,1), xlim=(0,maximum(plot_time)),
                yticks=(0:0.1:1, 0:0.1:1))
        end
        savefig(string(pwd(), "\\results\\$(sim_numbers)_fractions_over_time.png"))
    end
end
