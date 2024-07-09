function plot_active_mass_fractions(sim_number, time_index)
    """
    This function makes a horizontal bar plot of the active mass fractions.

    Arguments:
    sim_number:         The simulation number
    time_index:         The index (in time) of the data that needs to be plot_concentration_domain
                        If time_index = -1, take the most recent data

    Returns:
    A barplot (and saves it)
    """

    #Get stored data
    if isfile(string(pwd(), @sprintf("\\results\\%04d\\results1D.jld2", sim_number)))
        directory = string(pwd(), @sprintf("\\results\\%04d\\results1D.jld2", sim_number))
    else
        directory = string(pwd(), @sprintf("\\results\\%04d\\results2D.jld2", sim_number))
    end

    bac_vecint, bac_matfloat, bac_matint, bac_matbool = load(directory, "bac_saved_vecint", "bac_saved_matfloat", "bac_saved_matint", "bac_saved_matbool")
    
    # Get most recent data
    if time_index == -1
        time_index = findlast(bac_vecint.nBacs .!= 0)
    end

    # Extract required data
    nbac = bac_vecint.nBacs[time_index]
    actives = bac_matbool.active[time_index,1:nbac]
    species = bac_matint.species[time_index,1:nbac]
    radius = bac_matfloat.radius[time_index,1:nbac]

    # Calculate number of actives and their mass (and fractions)
    actives_per_species = [sum(actives .& (species .== i)) for i in 1:maximum(species)]
    total_moles_per_specie = [sum((radius[actives .& (species .== i)]).^3 .* (4 ./ 3) .* pi) for i in 1:maximum(species)] .* 500000 ./ 22.6232 
    total_moles = sum(total_moles_per_specie)
    mole_fractions = total_moles_per_specie ./ total_moles

    # Names and colors for Anaerobic Fermentation
    bac_names = ["BO" "AM" "HM"]
    color_list = ["#A020F0" "#E8BB00"  "#3bd6c6"]

    groupedbar(transpose(mole_fractions),
                bar_position= :stack,
                orientation = :horizontal,
                label= bac_names,
                bar_width=0.2,
                yticks=nothing,
                ylim=(0.7,1.3),
                color = color_list)
    
    # Vertical lines to indicate true division over the species           
    vline!([1/3, 2/3],
            linestyle=:dash,
            linewidth = 2,
            linecolor = :black,
            label=nothing)

    savefig(string(pwd(), "\\results\\$(sim_number)\\mass_fractions at $((time_index-1) * 12) h .png"))
end
