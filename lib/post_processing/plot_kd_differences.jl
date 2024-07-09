function plot_kd_differences()
    """
    This function is used to plot the active mass of AM to investigate the change in kd from fixed to relative values
    """

    # Select directory
    directory1 = string(pwd(), @sprintf("\\results\\9932\\results2D.jld2"))
    directory2 = string(pwd(), @sprintf("\\results\\9920\\results2D.jld2"))
    both_direct = [directory1, directory2]

    # Create clean plot
    plot()

    # Names and labels
    bac_names = ["BO" "AM" "HM"]
    labels = ["a = fixed", "a = 10%"]
    linestyles = [:solid, :dot]

    for direc in 1:2
        # Get data
        bac_vecint, bac_matfloat, bac_matint, bac_matbool = load(both_direct[direc], "bac_saved_vecint", "bac_saved_matfloat", "bac_saved_matint", "bac_saved_matbool")
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
    
            # Caluclate moles and fractions for each specie
            total_moles_per_specie = [sum((radius[actives .& (species .== i)]).^3 .* (4 ./ 3) .* pi) for i in 1:maximum(species)] .* 500000 ./ 22.6232 
            total_moles = sum(total_moles_per_specie)
            mole_fractions = total_moles_per_specie ./ total_moles
            store_mass[:,j] = total_moles_per_specie
            store_fractions[:,j] = mole_fractions
        end

        plot_time = (0:lastnonzero-1) * 12

        # Plotting the data
        plot!(plot_time, store_mass[2,:], color = "#E8BB00", label=labels[direc], linestyle=linestyles[direc], linewidth=2,
        xlabel="Time (h)", ylabel="Mass (mole)", xticks=(0:24:maximum(plot_time), 0:24:maximum(plot_time)),
        xlim=(0,125), ylim=(0, 3*(10^-12)))

    end
    savefig(string(pwd(), "\\results\\kd_differences.png"))
end
