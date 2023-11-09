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
        # Detachment based on detachment time and timestep (and size)
        include(string(pwd(), "\\lib\\determine_where_bacteria_in_grid.jl"))
        include(string(pwd(), "\\lib\\detachment\\calcTimeOfDetach.jl"))
        include(string(pwd(), "\\lib\\bacteria\\killBacs.jl"))

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

        # If timestep higher than detachment time --> detachment
        detachment_index = findall(1 .<= ratio .< Inf)
        bac_detach = zeros(length(detachment_index) * size(grid2bac, 3)) # Space for when all those grid cells are maximal occupied
        nDetach = 0 # Number so far
        for k in eachindex(detachment_index)
            y_index = detachment_index[k][1]
            x_index = detachment_index[k][2]
            grid_position = findall(grid2bac[y_index, x_index, :] .!= 0) # Finds which indices in the third dimension are filled
            # iBacs selects the actual values in the third dimension, these values are the bacterial indices in the bac struct
            iBacs = grid2bac[y_index, x_index, grid_position]
            n_temp = length(iBacs) # Number found this round
            bac_detach[nDetach + 1 : nDetach + n_temp] = iBacs # Store the indices in corresponding location
        end

        bac_detach = bac_detach[1:nDetach] # remove non-used places
        nCellsDetach = length(bac_detach)

        if nCellsDetach > 0
            bac = killBacs!(bac, bac_detach)
        end

        # Due to erosion, some bacteria on the outside are way too small 
        # Only factor 2 smaller than inactive bacteria should be reached
        # cells within the granule will not be removed
        mask_tooSmall = bac.MolarMass * constants.bac_MW .< constants.min_bac_mass_grams / 2
        xcentre = mean(bac.x[bac.active])
        ycentre = mean(bac.y[bac.active])
        # Why not this like before??
        # xcentre = grid.dx * grid.nx / 2
        # ycentre = grid.dy * grid.ny / 2
        dist = sqrt.((bac.x .- xcentre) .^2 + (bac.y .- ycentre) .^2)

        # A cell is considered on the "outside" when the difference in distance to the centre
        # compared to the distance of the furthest active and large enough cell is larger than 
        # 1 boundary layer thickness
        mask_outside = dist .> max(dist[bac.active .& .!mask_tooSmall]) - grid.blayer_thickness
        mask_outsideCellRemoval = mask_tooSmall .& mask_outside
        nCellsRemoved = sum(mask_outsideCellRemoval)

        if nCellsRemoved > 0
            bac = killBacs!(bac, mask_outsideCellRemoval)
            nCellsDetach = nCellsDetach + nCellsRemoved
        end


    # elseif settings.detachment == "suspension"
    elseif settings.model_type == "suspension"
        # Remove cells when they grow slower than dilution rate
        mask_remove = bac.mu .< (1 / invHRT)
        nCellsDetach = sum(mask_remove)

        if nCellsDetach > 0
            bac = killBacs!(bac, mask_remove)
        end
    else
        throw(ErrorException("Detachment of $(settings.detachment) is unknown"))
    end

    println("$(nCellsDetach) cells removed from the granule")
    return bac
end

