function bacteria_detachment!(bac, grid, constants, settings, timestep, invHRT)
    """
    This function implements rough detachment: remove all bacteria that are outside of maximum radius of granule
    bac is a struct containing all information regarding the bacteria
    constants is a struct containing all simulation constants
    settings is a struct contains the simulation parameters
    """
    println(pwd())

    if settings.detachment in ("none", "SBR")
        return

    elseif settings.detachment == "naive"
        include(string(pwd(), "\\lib\\bacteria\\killBacs.jl"))
        # Detach based on distance from centre
        bac_distance_from_centre = sqrt.((bac.x .- grid.dx * grid.nx / 2) .^2 + (bac.y .- grid.dy * grid.ny / 2) .^2)
        bac_detach = bac_distance_from_centre .> constants.max_granule_radius # Booleans
        nCellsDetach = sum(bac_detach)

        if nCellsDetach > 0
            bac = killBacs!(bac, bac_detach)
        end

    elseif settings.detachment == "mechanistic"
        include(string(pwd(), "\\lib\\determine_where_bacteria_in_grid.jl"))
        include(string(pwd(), "\\lib\\detachment\\calcTimeOfDetach.jl"))

        # grid2bac is a (nx * ny * ?) matrix containing the bacteria in each grid cell
        # grid2nBacs is a (nx * ny) matrix containing the number of bacteria in each grid cell
        # An update of these matrices is required as a division has occurred 
        grid2bac, grid2nBacs = determine_where_bacteria_in_grid(grid, bac)
        T = calcTimeOfDetach(bac, grid, grid2bac, grid2nBacs, constants)

        # Determine ratio of timestep of detachment time
        ratio = timestep ./ T

        # If timestep lower than detachment time --> decrease mass (erosion) (radius is updated when needed)
        erosion_index = findall(0 .< ratio .< 1)
        for k in eachindex(erosion_index)
            y_index = erosion_index[k][1]
            x_index = erosion_index[k][2]
            r = ratio[erosion_index[k]]
            grid_position = findall(grid2bac[y_index, x_index, :] .!= 0) # Finds which indices in the third dimension are filled
            # iBacs selects the actual values in the third dimension, these values are the bacterial indices in the bac struct
            iBacs = grid2bac[y_index, x_index, grid_position]

            for bac_index in iBacs
                # Mass gets adjusted according to the time that has passed and the detachment Time
                bac.molarMass[bac_index] = bac.molarMass[bac_index] * (1 - r)
            end
        end

        # Continue here

    end

end

