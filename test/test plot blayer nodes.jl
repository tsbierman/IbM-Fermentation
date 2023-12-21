# Getting what is required
using Plots

function circleShape(h, k, r)
    """
    Function to plot circles
    """
    thet = LinRange(0, 2*pi, 500)
    return h .+ r * sin.(thet), k .+ r*cos.(thet)
end

# Find files
include(string(pwd(), "\\lib\\pre_processing\\create_mat.jl"))
include(string(pwd(), "\\lib\\determine_diffusion_region.jl"))
include(string(pwd(), "\\lib\\determine_where_bacteria_in_grid.jl"))

# Read_file
filename = string(pwd(), "\\test\\test_file.xlsx")

# Initilising 
grid, bac, constants, settings, init_params = create_mat(filename)
grid2bac, grid2nBacs = determine_where_bacteria_in_grid(grid, bac)

# Plot all grid points
nodeEndCoordinatesX = (1:grid.nx) * grid.dx
nodeEndCoordinatesY = (1:grid.ny) * grid.dy
# Makes Centered
nodeMidCoordinatesX = nodeEndCoordinatesX .- grid.dx/2
nodeMidCoordinatesY = nodeEndCoordinatesY .- grid.dy/2
# Make all coordination possible
plot_centre_X = repeat(nodeMidCoordinatesX, length(nodeMidCoordinatesY))
plot_centre_Y = repeat(nodeMidCoordinatesY, inner=length(nodeMidCoordinatesX))

# Get nodes that are estimated to be diffusion nodes
x_Nodes, y_Nodes = getDiffusionNodes(bac, grid)
x_init_diff = nodeMidCoordinatesX[x_Nodes]
y_init_diff = nodeMidCoordinatesY[y_Nodes]

# All combinations of those
plot_init_diff_X = repeat(x_init_diff, length(y_init_diff))
plot_init_diff_Y = repeat(y_init_diff, inner=length(x_init_diff))

# Offset to check whether initial diffusion nodes are within reach
offsetX = grid.blayer_thickness + grid.dx
offsetY = grid.blayer_thickness + grid.dy

x_min = minimum(bac.x) - offsetX
x_max = maximum(bac.x) + offsetX
y_min = minimum(bac.y) - offsetY
y_max = maximum(bac.y) + offsetY

# Get actual diffusion nodes
diffusion_region, focus_region = determine_diffusion_region(grid2bac, grid2nBacs, bac, grid)
trans_diff = transpose(diffusion_region)
flattened = Base.Iterators.flatten(trans_diff)
# Select correct nodes
plot_diff_X = plot_centre_X[collect(flattened)]
plot_diff_Y = plot_centre_Y[collect(flattened)]

# Plotting for both granule and suspension mode
if settings.model_type == "granule"
    scatter(plot_centre_X, plot_centre_Y, mc=:black, ms = 0.6, ylimits = (2.300e-4,2.85e-4), xlimits = (2.300e-4,2.85e-4)) # plot grid centres (limits for granule)
    for index in eachindex(bac.x)
        plot!(circleShape(bac.x[index], bac.y[index], bac.radius[index]), seriestype = [:shape,],
        c=:green, fillalpha=0.5, linealpha=0.5, legend = false, aspect_ratio=1)
    end
    scatter!(plot_init_diff_X, plot_init_diff_Y, ms=2, mc=:red, ylimits = (2.300e-4,2.85e-4), xlimits = (2.300e-4,2.85e-4)) # plot initial diffusion nodes
    # Test whether they are all within offset
    hline!([y_min, y_max])
    vline!([x_min, x_max])
    scatter!(plot_diff_X, plot_diff_Y, mc=:blue, ms = 2, ylimits = (2.300e-4,2.85e-4), xlimits = (2.300e-4,2.85e-4)) # plot grid centres (limits for granule)
    
    for index in eachindex(bac.x)
        plot!(circleShape(bac.x[index], bac.y[index], grid.blayer_thickness), seriestype = [:shape,], 
        legend = false, aspect_ratio=1, fillalpha=0.1, linealpha=0.1)
    end
    # savefig("Granule_diffusion_nodes.png")

elseif settings.model_type == "suspension"
    scatter(plot_centre_X, plot_centre_Y, mc=:black, ms = 0.5, fillalpha=0.1) # plot grid centres 
    for index in eachindex(bac.x)
        plot!(circleShape(bac.x[index], bac.y[index], bac.radius[index]), seriestype = [:shape,],
        c=:green, fillalpha=0.5, linealpha=0.5, legend = false, aspect_ratio=1)
    end
    scatter!(plot_init_diff_X, plot_init_diff_Y, ms=0.2, mc=:red, fillalpha=0.1,
    ylimits= (y_min - 0.1e-4, y_max + 0.1e-4), xlimits = (x_min - 0.1e-4, x_max + 0.1e-4))
    # Test whether they are all within offset
    hline!([y_min, y_max])
    vline!([x_min, x_max])
    scatter!(plot_diff_X, plot_diff_Y, mc=:blue, ms = 2)

    for index in eachindex(bac.x)
        plot!(circleShape(bac.x[index], bac.y[index], grid.blayer_thickness), seriestype = [:shape,], 
        legend = false, aspect_ratio=1, fillalpha=0.1, linealpha=0.1)
    end
    # savefig("Suspension_diffusion_nodes_all.png")
end
plot!(aspect_ratio = 1, legend=false)
