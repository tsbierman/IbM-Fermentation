function rand_circle(N, x_centre, y_centre, r)
""" 
This function takes the amount of candidates required and generates a set
of x and y coordinates within the circle.
"""
    # Create estimate that is actually in the circle
    Ns = round(4 / pi * N + 2.5 * sqrt(N) + 100)
    X = rand(1, Ns) * (2*r) .- r # Random decimal time Diameter, shifted to be around 0.
    Y = rand(1, Ns) * (2*r) .- r
    I = findall(sqrt.(X .^2 .+ Y .^2) <= r) # Check which are within radius
    X = X[I[1:N]] .+ x_centre # Select and move to existing centre
    Y = Y[I[1:N]] .+ y_centre

    # Alternative: This generates only points that are already in the circle, no need to check the distance.
    # radius = r * sqrt(rand(N,1))
    # theta = 2 * pi * rand(N,1)
    # X = radius .* cos(theta) + x_centre
    # Y = radius .* sin(theta) + y_centre
    return X, Y
end


function blue_noise_circle(n, x_centre, y_centre, r)
"""
This function evaluates coordinate candidates based on distance to the existing points.
It stores the set that is the furthest away. This is done for n sets of coordinates.
"""
    m = 20 # Amount of potential candidates that will be generated
    X = zeros(n, 1)
    Y = zeros(n, 1)
    X[1] = x_centre
    Y[1] = y_centre
    npoints = 1

    for cell in 1:(n-1)
        # Generate candidate coordinates
        candidate_x, candidate_y = rand_circle(m, x_centre, y_centre, r)
        # Distance all candidates and existing points
        distance = sqrt.((X[1:npoints] .- candidate_x) .^2 .+ (Y[1:npoints] .- candidate_y) .^2) 
        # First find closest existing point for each candidate, then index of furthest existing point over all candidates.
        _, index = findmax(minimum(distance, dims = 1))

        # Store
        npoints = npoints + 1
        X[npoints] = candidate_x[index]
        Y[npoints] = candidate_y[index]
    end

    return X, Y
end


function distribute_microcolonies(nColonies, nBacPerCol, r_colony, xrange, yrange)
    nsections = ceil(sqrt(nColonies) * 1.1)
    xlist = range(xrange[1], xrange[2], nsections)
    ylist = range(yrange[1], yrange[2], nsections)

    space_margin = xlist[2] - xlist[1]


end

#Initialise structs that will have to be used later
# Structs need to be declared at top level

loading_file = string(Base.source_dir(), "\\","loadPresetFile.jl")
include(loading_file)

struct General
    properties::Dict{Symbol, Any}
end
General() = General(Dict{Symbol, Any}())

Base.getproperty(x::General, property::Symbol) = getfield(x, :properties)[property]
Base.setproperty!(x::General, property::Symbol, value) = getfield(x, :properties)[property] = value
Base.propertynames(x::General) = keys(getfield(x, :properties))

import XLSX

# This needs to be changed if the test_file is moved to another directory.
code_folder = dirname(dirname(Base.source_dir()))
filename = string(code_folder, "\\planning\\test_file.xlsx")

println(">>>>>>>>>>>>>>>>>> LOADING EXCEL FILE")
grid, bac_init, constants, settings, init_params = loadPresetFile(filename)

# Molar mass is 60% of maximum weigth
molarMass = 0.6 * constants.max_bac_mass_grams / constants.bac_MW
radius = ( (molarMass * constants.bac_MW / constants.bac_rho) * (3 / (4 * pi))) ^ (1/3)

println(">>>>>>>>>>>>>>>> INITIALISING BACTERIA")

bac = General()

if settings.model_type in ("granule", "mature granule")

    bac.x, bac.y = blue_noise_circle(bac_init.start_nBac, grid.nx / 2 * grid.dx, grid.ny / 2 * grid.dy, bac_init.granule_radius)

elseif settings.model_type in ("suspension")

    margin = 0.2 * grid.dx * grid.nx # 20% of simulation domain as margin for letting suspensions growth (empirical)
    xrange = [margin, grid.dx*grid.nx - margin]
    yrange = xrange # assume square domain

    r_colony = (bac_init.start_nBacPerColony * radius * constants.kDist) / 5 # Empirical, 1/5 * diameter if all cell next to each other.
    bac.x, bac.y = distribute_microcolonies(bac_init.start_nColonies, bac_init.start_nBacPerColony, r_colony, xrange, yrange)

end
