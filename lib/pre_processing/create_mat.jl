using Printf

function rand_circle(N, x_centre, y_centre, r)
    """ 
    This function generates a set of coordinates around a centre, at most r away.

    Arguments
    N:                      Amount of coordinate sets that have to be generated
    x_centre, y_centre      X and Y coordinate around which the colony has to be generated
    r:                      Maximum radius of the colony

    Returns
    X, Y                    Vectors [N,] containing X and Y coordinates of the candidates
    """
    # Create estimate that is actually in the circle
    Ns = round(4 / pi * N + 2.5 * sqrt(N) + 100)
    X = rand(1, Int(Ns)) * (2*r) .- r           # Random decimal times Diameter, shifted to be around 0.
    Y = rand(1, Int(Ns)) * (2*r) .- r           # Random decimal times Diameter, shifted to be around 0.
    
    I = findall(sqrt.(X .^2 .+ Y .^2) .<= r)    # Check which are within radius
    X = X[I[1:N],:] .+ x_centre                 # Select and move to existing centre
    Y = Y[I[1:N],:] .+ y_centre                 # Select and move to existing centre

    X = reshape(X, 1, :)
    Y = reshape(Y, 1, :)

    # Alternative: This generates only points that are already in the circle, no need to check the distance.
    # radius = r * sqrt(rand(N,1))
    # theta = 2 * pi * rand(N,1)
    # X = radius .* cos(theta) + x_centre
    # Y = radius .* sin(theta) + y_centre
    return X, Y
end


function blue_noise_circle(n, x_centre, y_centre, r)
    """
    This function calls rand_circle to generate candidate coordinates.
    It then evaluates this coordinates and chooses to add the one with the 
    maximum minimum distance, so furthest away from any of the exisiting points.

    Arguments
    n:                      The amount of points (bacteria) that need to be generated
    x_centre, y_centre      The x and y centre around which points need to be generated
    r:                      Maximum distance of points from the centre.

    Returns
    X, Y                    Vectors [n,] containing the coordinates of the bacteria
    """
    m = 20 # Amount of potential candidates that will be generated
    X = zeros(n)
    Y = zeros(n)
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
    """
    This function creates several microcolonies over the available grid. 
    It does so by calling the blue_noise_circle function to create a cluster 
    of microbial cells at every starting point.

    Arguments
    nColonies:          The amount of colonies to be generated
    nBacPerCol:         The amount of bacteria per colony
    r_colony:           The allowed radius of a colony
    xrange:             Minimum and maximum x-coordinate allowed
    yrange:             Minimum and maximum y-coordinate allowed

    Returns
    x, y                Vectors [nColonies*nBacPerCol] containing all x and y coordinates, 
                        respectively, of the bacteria
    """
    # Estimate how many colonies per axis and create equally spaced locations for the colonies
    nsections = Int(ceil(sqrt(nColonies) * 1.1))
    xlist = transpose(range(xrange[1], xrange[2], nsections))
    ylist = transpose(range(yrange[1], yrange[2], nsections))

    # Calculate distance between 2 colonies
    space_margin = xlist[2] - xlist[1]

    # Generate 2 vectors that combined include every possible set of coordinates
    x = transpose(repeat(xlist, 1, nsections))
    y = transpose(repeat(ylist, inner=(1,nsections)))

    # Create some noise in both coordinates and update the vectors
    noise_x = rand(size(x,1), size(x,2)) * 0.6 * space_margin .- 0.3 * space_margin
    noise_y = rand(size(y,1), size(y,2)) * 0.6 * space_margin .- 0.3 * space_margin

    x = x + noise_x
    y = y + noise_y

    # Shuffle the coordinates and only select how many we will need
    i = randperm(length(x))
    i = i[1:nColonies]

    x = x[i]
    y = y[i]
    
    # For every coordinate (1 cell), generate a colony around the cell and store its coordinates    
    for index in eachindex(x)
        x_microcol, y_microcol = blue_noise_circle(nBacPerCol, x[index], y[index], r_colony)
        x = [x;x_microcol[2:end]]
        y = [y;y_microcol[2:end]]
    end
    return x, y
end


function AMXinside(bac, grid, constants)
    """
    This function assigns all the bacteria a specie. The closest to the centre are assigned with AMX,
    the others are random. The amount assigned with AMX remains random.
    This function is called in case of a mature granule and is Nitrospira specific.

    Arguments
    bac:                A "General" struct containing all parameters related to the bacteria
    grid:               A "General" struct containing all parameters related to the grid
    constants:          A "General" struct containing all the simulation constants

    Returns
    species:            A vector (nBac,) containing a number, representing the species of the organism
    """
    # Calculate distance from the centre and sort based on distance
    distance = sqrt.((bac.x .- (grid.nx / 2 * grid.dx)) .^2 + (bac.y .- (grid.ny / 2 * grid.dy)) .^2)
    I = sortperm(distance)   # Returns indices

    # Determine index of AMX specie. Randomly assign species to bacteria and see how many are AMX
    iAMX = findall(constants.speciesNames .== "B2")[1]
    nAMX = sum(rand(1:length(constants.speciesNames), size(bac.x)) .== iAMX)

    # Generate storage for species and select nAMX closest to centre to be AMX
    species = zeros(size(bac.x))
    species[I[1:nAMX]] .= iAMX

    # Generate random specie index for all the bacteria that do not have a specie (so 0 in species)
    # 1 less so we miss a specie. Then everything equal or higher than iAMX will shift upwards
    species_other = rand(1:length(constants.speciesNames) .- 1, size(species[species .== 0]))
    species_other[species_other .>= iAMX] = species_other[species_other .>= iAMX] .+ 1
    species[I[nAMX+1:end]] = species_other
    return species
end

function shoving_loop(bac, grid, constants, n)
    """
    This function calls the bacteria_shove multiple thermodynamic_parameters
    
    Arguments
    bac:                A "General" struct containing all parameters related to the bacteria
    grid:               A "General" struct containing all parameters related to the grid
    constants:          A "General" struct containing all the simulation constants
    n:                  Amount of times the shoving algorithm is called

    Returns
    bac:                Bac struct with updates x and y coordinates
    """
    for gg in 1:n
        bac = bacteria_shove!(bac, grid, constants)
    end
    return bac
end

function create_mat(filename, simulation_number)
    """
    This function reads an excel file and extracts all pre-set parameters from it.
    It stores these parameters in structs, which can be used later in the simulation.
    Simulation_number guides how this function is used. The normal use case is a call
    with a simulation number from 1 to 9999. This will create a .jld2 file with the variables saved.
    The user has to call IbM(simulation_number) themselves. 
    When 0 is used as simulation_number, create_mat will
    still save everything, but will call IbM(0) directly. This use case is for starting up and
    should only be called with the file (start_up.xlsx).
    When a negative number is supplied, this function will not save the variables but will
    return them. This is useful when testing as saving might not be desired then.

    Arguments:
    filename:           An excel file containing all parameters
    simulation_number   The number the user wants this simulation stored under.

    Returns
    grid:               A "General" struct containing all parameters related to the grid
    bac:                A "General" struct containing all parameters related to the bacteria
    constants:          A "General" struct containing all the simulation constants
    settings:           A "General" struct containing all the settings of the simulation
    init_params:        A "General" struct containing the parameters values at the start of the simulation
    """

    println(">>>>>>>>>>>>>>>>>> LOADING EXCEL FILE")

    grid, bac_init, constants, settings, init_params = loadPresetFile(filename)

    # Initial Molar Mass is 60% of maximum molar Mass
    molarMass = 0.6 * constants.max_bac_mass_grams / constants.bac_MW                           # [mol]
    radius = ( (molarMass * constants.bac_MW / constants.bac_rho) * (3 / (4 * pi))) ^ (1/3)     # [m]

    println(">>>>>>>>>>>>>>>> INITIALISING BACTERIA")

    bac = General()

    if settings.model_type in ("granule", "mature granule")

        # Create a single colony of bacteria around the centre
        bac.x, bac.y = blue_noise_circle(bac_init.start_nBac, grid.nx / 2 * grid.dx, grid.ny / 2 * grid.dy, bac_init.granule_radius)

    elseif settings.model_type in ("suspension",)

        margin = 0.2 * grid.dx * grid.nx                # 20% of simulation domain as margin for letting suspensions growth (empirical)
        xrange = [margin, grid.dx*grid.nx - margin]
        yrange = xrange                                 # assume square domain

        # Create several colonies with some bacteria each
        r_colony = (bac_init.start_nBacPerColony * radius * constants.kDist) / 5 # Empirical, 1/10 * diameter if all cell next to each other.
        bac.x, bac.y = distribute_microcolonies(bac_init.start_nColonies, bac_init.start_nBacPerColony, r_colony, xrange, yrange) # Generate all coordinates
    end

    # Set parameters for every of the bacteria
    bac.molarMass = ones(length(bac.x)) * molarMass         # [mol]
    bac.radius = ones(length(bac.x)) * radius               # [m]
    bac.active = BitArray(ones(size(bac.x)))                # Binary/Boolean

    # Shove bacteria to prevent overlapping at the start. The 5 is arbritrary.
    bac = shoving_loop(bac, grid, constants, 5)

    if settings.model_type in ("granule", "mature granule")
        # Remove bacteria that are outside the maximum granule radius due to shoving
        keep = sqrt.((bac.x .- (grid.dx * grid.nx / 2)) .^2 + (bac.y .- (grid.dy * grid.ny / 2)) .^2 ) .<= bac_init.granule_radius
        println("$(size(bac.x, 1)- sum(keep)) Bacteria removed outside of starting granule")
        bac.x = bac.x[keep]
        bac.y = bac.y[keep]
        bac.radius = bac.radius[keep]
        bac.molarMass = bac.molarMass[keep]
        bac.active = bac.active[keep]
    end

    if settings.model_type == "mature granule"
        bac.species = AMXinside(bac, grid, constants) # If granule already mature, assign AMX types on the inside (Nitrospira specific)
    else
        bac.species = rand((1:length(constants.speciesNames)), size(bac.x)) # Random species
    end

    println("$(length(bac.x)) starting bacteria in the system")

    for specie in eachindex(constants.speciesNames)
        println("\t $(sum(bac.species .== specie)) $(constants.speciesNames[specie])")
    end

    if simulation_number < 0 
        # For testing purposes
        println(">>>>>>>>>>>>>> DONE LOADING!")
        return grid, bac, constants, settings, init_params

    elseif simulation_number == 0 
        # Run activation simulation
        results_file = @sprintf("sim_%04d.jld2", simulation_number)
        # Save
        save(results_file, "grid", grid, "bac", bac, "constants", constants, "settings", settings, "init_params", init_params)
        println(">>>>>>>>>>>>>> DONE LOADING AND SAVING!")
        # Directly call IbM
        IbM(simulation_number)

    else 
        # Normal use case
        results_file = @sprintf("sim_%04d.jld2", simulation_number)
        save(results_file, "grid", grid, "bac", bac, "constants", constants, "settings", settings, "init_params", init_params)
        println(">>>>>>>>>>>>>> DONE LOADING AND SAVING!")

    end
end

# Import everything necessary while still in top-level scope
import XLSX
using Random
using Plots
using InvertedIndices
using Statistics
using DSP
using DifferentialEquations
using ODE
using DataStructures
using FileIO
using TickTock

# # Initialise struct "General" that can dynamically add properties
# include(string(pwd(), "\\lib\\Struct_Module.jl"))

# # Include required files
# include(string(pwd(), "\\lib\\pre_processing\\loadPresetFile.jl"))
# include(string(pwd(), "\\lib\\bacteria\\bacteria_shove.jl"))
# include(string(pwd(), "\\IbM.jl"));
