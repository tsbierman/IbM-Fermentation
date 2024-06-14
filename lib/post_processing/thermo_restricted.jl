function circleShape(h, k, r)
    """
    Function to plot circles
    """
    thet = LinRange(0, 2*pi, 500)
    return h .+ r * sin.(thet), k .+ r*cos.(thet)
end

function HEX2RGB(c)

    rC, gC, bC = zeros(length(c)), zeros(length(c)), zeros(length(c))
    for idx in eachindex(c)
        color = c[idx]
        RGB = parse(UInt, color[2:end], base=16)
        rC[idx] = (RGB >> 16) & 0xFF
        gC[idx] = (RGB >> 8) & 0xFF
        bC[idx] = RGB & 0xFF
    end
    return rC, gC, bC
end

function thermo_restricted(sim_number, finished)
    """
    Shows which gridcells do not follow thermodynamic restrictions but do reaction nontheless
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
    bac_saved_vecint, bac_saved_matint, bac_saved_matfloat, conc_saved = load(directory, "bac_saved_vecint", "bac_saved_matint", "bac_saved_matfloat", "conc_saved")

    # Get most recent index, nBacs and all species
    lastnonzero = findlast(bac_saved_vecint.nBacs .!= 0)
    nBacs = bac_saved_vecint.nBacs[lastnonzero]
    species = bac_saved_matint.species[lastnonzero, 1:nBacs]

    short_bac_vecfloat = VectorFloat_struct()
    bac_x = bac_saved_matfloat.x[lastnonzero, 1:nBacs]
    bac_y = bac_saved_matfloat.y[lastnonzero, 1:nBacs]
    mu_s = bac_saved_matfloat.mu[lastnonzero, 1:nBacs]
    short_bac_vecfloat.x = bac_x
    short_bac_vecfloat.y = bac_y

    concentrations = conc_saved[lastnonzero, :, :, :]
    grid_float.blayer_thickness = grid_float.dx

    grid2bac, grid2nBacs = determine_where_bacteria_in_grid(grid_float, grid_int, short_bac_vecfloat)
    diffusion_region, focus_region = determine_diffusion_region(grid2bac, grid2nBacs, short_bac_vecfloat, grid_float, grid_int)

    # ORDER: (NH3 is not used in catabolisms)
    # Butyrate, Acetate, CH4, HCO3, H2, H2O, H

    Gibbs_database =    [-352.6, -369.4, -50.8, -586.9, 0, -237.2, 0]
    Enthalpy_database = [-535.0, -485.8, -74.8, -692.0, 0, -285.8, 0]
    Ts = 28 + 273.15        # K
    T = 35 + 273.15         # K
    R = 8.314e-3            # kJ/K/mol

    # Butyrate, Acetate, CH4, HCO3, H2, H2O, H
    BO_cat = [-1.0, 2.0, 0.0, 0.0, 2.0, -2.0, 1.0]
    AM_cat = [0.0, -1.0, 1.0, 1.0, 0.0, -1.0, 0.0]
    HM_cat = [0.0, 0.0, 0.25, -0.25, -1.0, 0.75, -0.25]

    # Collect all catabolisms
    all_cat = [BO_cat, AM_cat, HM_cat]

    # Calculate dG0
    BO_G0 = sum(BO_cat .* Gibbs_database)
    AM_G0 = sum(AM_cat .* Gibbs_database)
    HM_G0 = sum(HM_cat .* Gibbs_database)

    # Calculate dH0
    BO_H0 = sum(BO_cat .* Enthalpy_database)
    AM_H0 = sum(AM_cat .* Enthalpy_database)
    HM_H0 = sum(HM_cat .* Enthalpy_database)

    # Correct G0 for Temperature
    BO_G0T = BO_G0 * T/Ts + BO_H0 * (Ts-T)/Ts
    AM_G0T = AM_G0 * T/Ts + AM_H0 * (Ts-T)/Ts
    HM_G0T = HM_G0 * T/Ts + HM_H0 * (Ts-T)/Ts

    # Combine all the Corrected G0s
    G0Ts = [BO_G0T, AM_G0T, HM_G0T]

    # Create storage, 1 dimension for every Metabolism, plus another for the combination
    thermo_checks = zeros(size(grid2nBacs, 1), size(grid2nBacs, 1), 4)
    thermo_checks[diffusion_region, :] .= 2

    # Traverse over all grid cells
    for ix in axes(grid2nBacs, 2)
        for iy in axes(grid2nBacs, 1)

            # Only continue checking when the grid cell contains bacteria
            if grid2nBacs[iy, ix] > 0
                number_of_bacs = grid2nBacs[iy, ix]
                current_concs = [concentrations[iy,ix, 1:5]; 1; 1e-7] # Get concs of Butyrate, Acetate, CH4, CO2 and H2, added H2O to have no effect and pH 7 for H+

                # For every bacteria on that location
                for bac_num in 1:number_of_bacs
                    # Get bac index in whole system, its specie number and the catabolism
                    bac_index = grid2bac[iy,ix,bac_num]
                    bac_specie = species[bac_index]
                    used_cat = all_cat[bac_specie]

                    # Elevate all concentrations to their coefficient
                    adjusted_concs = prod(current_concs .^ used_cat)
                    # Correct the G0T for the real concentrations
                    G1 = G0Ts[bac_specie] + R * T * log(adjusted_concs)

                    # If G1 is positive, but active growth has been observed --> place a 1 in the matrix for that specie
                    if G1 >= 0 && mu_s[bac_index] > 0
                        thermo_checks[iy, ix, bac_specie] = 1
                        # println(current_concs)
                    # If G1 is negative OR positive, but with negative growth, this specie adheres to thermodynamics, place a -1
                    elseif G1 < 0 || (G1 >=0 && mu_s[bac_index] < 0)
                        thermo_checks[iy, ix, bac_specie] = -1
                    end
                end

                # If any specie is doing impossible on this grid cell, add a 1 to the total layer as well
                if any(thermo_checks[iy, ix, 1:3] .== 1)
                    thermo_checks[iy, ix, 4] = 1
                elseif any(thermo_checks[iy, ix, 1:3] .== -1)
                    thermo_checks[iy, ix, 4] =-1
                end
            end
        end
    end

    names = ["BO", "AM", "HM", "Total"]

    red_r, red_g, red_b = HEX2RGB(["#E74F4E"])
    green_r, green_g, green_b = HEX2RGB(["#00D400"])
    my_red = RGB(red_r[1]/255, red_g[1]/255, red_b[1]/255)
    my_green = RGB(green_r[1]/255, green_g[1]/255, green_b[1]/255)
    
    for name_index in eachindex(names)
        plot(aspect_ratio = 1)
        for ix in axes(grid2nBacs, 2)
            for iy in axes(grid2nBacs, 1)
                # x_coor = (ix * grid_float.dx - 0.5 * grid_float.dx) * 1e6
                # y_coor = (iy * grid_float.dy - 0.5 * grid_float.dy) * 1e6

                # if thermo_checks[iy, ix, name_index] .== 0 && grid2nBacs[iy, ix] > 0
                #     plot!(circleShape(x_coor, y_coor, 0.5 * grid_float.dx * 1e6), seriestype =[:shape], aspect_rate=1, c = my_green, linewidth = 0.1, legend = false, 
                #                       xlabel="Position along x-axis [µm]", ylabel="Position along y-axis [µm]")
                # elseif thermo_checks[iy, ix, name_index] .== 1
                #     plot!(circleShape(x_coor, y_coor, 0.5 * grid_float.dx * 1e6), seriestype =[:shape], aspect_rate=1, c = my_red, linewidth = 0.1, legend = false,
                #                       xlabel="Position along x-axis [µm]", ylabel="Position along y-axis [µm]")

                x_coor_left = ((ix-1) * grid_float.dx) * 1e6
                x_coor_right = (ix * grid_float.dx) * 1e6
                y_coor_bottom = ((iy-1) * grid_float.dy) * 1e6
                y_coor_up = (iy * grid_float.dy) * 1e6
                
                if thermo_checks[iy, ix, name_index] .== -1
                    plot!(Shape([(x_coor_right, y_coor_up), (x_coor_right, y_coor_bottom), (x_coor_left, y_coor_bottom), (x_coor_left, y_coor_up)]),
                    seriestype =[:shape], aspect_rate=1, c = my_green, linewidth = 0.1, legend = false, 
                                      xlabel="Position along x-axis [µm]", ylabel="Position along y-axis [µm]")
                elseif thermo_checks[iy, ix, name_index] .== 1
                    plot!(Shape([(x_coor_right, y_coor_up), (x_coor_right, y_coor_bottom), (x_coor_left, y_coor_bottom), (x_coor_left, y_coor_up)]), 
                    seriestype =[:shape], aspect_rate=1, c = my_red, linewidth = 0.1, legend = false,
                                      xlabel="Position along x-axis [µm]", ylabel="Position along y-axis [µm]")
                elseif thermo_checks[iy, ix, name_index] .== 2
                    plot!(Shape([(x_coor_right, y_coor_up), (x_coor_right, y_coor_bottom), (x_coor_left, y_coor_bottom), (x_coor_left, y_coor_up)]),
                    seriestype =[:shape], aspect_rate=1, c = :grey, linewidth = 0.1, legend = false, 
                                      xlabel="Position along x-axis [µm]", ylabel="Position along y-axis [µm]")
                end
            end
        end
        save_loc = string(pwd(), @sprintf("\\results\\%04d", sim_number))
        savefig("$(save_loc)\\$(names[name_index]).png")
    end
end
