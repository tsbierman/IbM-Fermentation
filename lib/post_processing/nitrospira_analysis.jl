function nitrospira_analysis(sim_number)
    """Does the analysis for the nitrospira validation cases."""
    
    # if sim_number in [9957, 9958, 9959]
    if sim_number in [9951, 9952, 9953]
        destination = string(pwd(), "\\results\\1.0_O2_mass_analysis.xlsx")
    elseif sim_number in [9954, 9955, 9956]
        destination = string(pwd(), "\\results\\1.5_O2_mass_analysis.xlsx")
    else
        destination = string(pwd(), "\\results\\3.0_O2_mass_analysis.xlsx")
    end

    # sim_row_dict = Dict(9957=>2,
    #                     9958=>3,
    #                     9959=>4,
    #                     9960=>2,
    #                     9961=>3,
    #                     9962=>4)

    sim_row_dict = Dict(9948=>2,
                        9949=>3,
                        9950=>4,
                        9951=>2,
                        9952=>3,
                        9953=>4,
                        9954=>2,
                        9955=>3,
                        9956=>4)

    bac_vecint, bac_matfloat, bac_matint, bac_matbool = load(string(pwd(), "\\results\\$(sim_number)\\results1D.jld2"), "bac_saved_vecint", "bac_saved_matfloat", "bac_saved_matint", "bac_saved_matbool")
    final_nbac = bac_vecint.nBacs[end-1]
    actives = bac_matbool.active[end-1,1:final_nbac]
    species = bac_matint.species[end-1,1:final_nbac]
    radius = bac_matfloat.radius[end-1,1:final_nbac]
    actives_per_species = [sum(actives .& (species .== i)) for i in 1:maximum(species)]
    total_moles_per_specie = [sum((radius[actives .& (species .== i)]).^3 .* (4 ./ 3) .* pi) for i in 1:maximum(species)] .* 500000 ./ 24.6 

    XLSX.openxlsx(destination, mode="rw") do xf
        sheet = xf[1]
        row = sim_row_dict[sim_number]
        sheet["B$(row)"] = sum(actives)
        sheet["C$(row)"] = actives_per_species

        sheet2 = xf[2]
        sheet2["B$(row)"] = total_moles_per_specie
    end
end
