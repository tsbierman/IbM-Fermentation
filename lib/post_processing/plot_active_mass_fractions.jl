function plot_active_mass_fractions(sim_number, time_index)

    if isfile(string(pwd(), @sprintf("\\results\\%04d\\results1D.jld2", sim_number)))
        directory = string(pwd(), @sprintf("\\results\\%04d\\results1D.jld2", sim_number))
    else
        directory = string(pwd(), @sprintf("\\results\\%04d\\results2D.jld2", sim_number))
    end

    bac_vecint, bac_matfloat, bac_matint, bac_matbool = load(directory, "bac_saved_vecint", "bac_saved_matfloat", "bac_saved_matint", "bac_saved_matbool")
    
    if time_index == -1
        time_index = findlast(bac_vecint.nBacs .!= 0)
    end

    nbac = bac_vecint.nBacs[time_index]
    actives = bac_matbool.active[time_index,1:nbac]
    species = bac_matint.species[time_index,1:nbac]
    radius = bac_matfloat.radius[time_index,1:nbac]
    actives_per_species = [sum(actives .& (species .== i)) for i in 1:maximum(species)]
    total_moles_per_specie = [sum((radius[actives .& (species .== i)]).^3 .* (4 ./ 3) .* pi) for i in 1:maximum(species)] .* 500000 ./ 22.6232 
    total_moles = sum(total_moles_per_specie)
    mole_fractions = total_moles_per_specie ./ total_moles

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
    
    vline!([1/3, 2/3],
            linestyle=:dash,
            linewidth = 2,
            linecolor = :black,
            label=nothing)

    savefig(string(pwd(), "\\results\\$(sim_number)\\mass_fractions at $((time_index-1) * 12) h .png"))
end
