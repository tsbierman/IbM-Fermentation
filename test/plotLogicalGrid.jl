function circleShape(h, k, r)
    """
    Function to plot circles
    """
    thet = LinRange(0, 2*pi, 500)
    return h .+ r * sin.(thet), k .+ r*cos.(thet)
end

function plotLogicalGrid(grid, logicalGrid)
    """
    Plots logical grid
    """
    # Plot all grid points
    nodeEndCoordinatesX = (1:grid.nx) * grid.dx
    nodeEndCoordinatesY = (1:grid.ny) * grid.dy
    # Makes Centered
    nodeMidCoordinatesX = nodeEndCoordinatesX .- grid.dx/2
    nodeMidCoordinatesY = nodeEndCoordinatesY .- grid.dy/2
    # Make all coordination possible
    plot_centre_X = repeat(nodeMidCoordinatesX, length(nodeMidCoordinatesY))
    plot_centre_Y = repeat(nodeMidCoordinatesY, inner=length(nodeMidCoordinatesX))

    trans_logical = transpose(logicalGrid)
    flattened = Base.Iterators.flatten(trans_logical)
    # Select correct nodes
    plot_diff_X = plot_centre_X[collect(flattened)]
    plot_diff_Y = plot_centre_Y[collect(flattened)]

    scatter(plot_centre_X, plot_centre_Y, mc=:black, ms = 0.6, ylimits = (2.300e-4,2.85e-4), xlimits = (2.300e-4,2.85e-4)) # plot grid centres (limits for granule)
    for index in eachindex(bac.x)
        plot!(circleShape(bac.x[index], bac.y[index], bac.radius[index]), seriestype = [:shape,],
        c=:green, fillalpha=0.5, linealpha=0.5, legend = false, aspect_ratio=1)
    end

    for index in eachindex(bac.x)
        plot!(circleShape(bac.x[index], bac.y[index], grid.blayer_thickness), seriestype = [:shape,], 
        legend = false, aspect_ratio=1, fillalpha=0.1, linealpha=0.1)
    end
    
    scatter!(plot_diff_X, plot_diff_Y, mc=:blue, ms = 2, ylimits = (2.300e-4,2.85e-4), xlimits = (2.300e-4,2.85e-4)) # plot grid centres (limits for granule)
    return plot!(legend=false)
end
