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


function muratio(mu, species, inc)
    if inc == "NoAlpha"
        return ones(size(mu))
    else
        println("Sorry, mu/mumax not implemented")
    end
end


function save_plot(i, xlim, ylim, bac_vecint, bac_matfloat, bac_matint, bac_matbool, bacNames, dT_save, sim_number, only_last, include_radius)

    nBacs = bac_vecint.nBacs[i]
    x = bac_matfloat.x[i, 1:nBacs] * 1e6            # In µm
    y = bac_matfloat.y[i, 1:nBacs] * 1e6            # In µm
    radius = bac_matfloat.radius[i, 1:nBacs] * 1e6  # In µm
    species = bac_matint.species[i, 1:nBacs]
    active = bac_matbool.active[i, 1:nBacs]
    mu = bac_matfloat.mu[i, 1:nBacs]

    # Calculus of mu/max(mu)
    # Recommended: 1.0 - 2.0  (inc = 'NoAlpha' -> alpha = 1; inc = 0 -> alpha = mu/max_mu)
    inc = "NoAlpha"
    muAlpha = muratio(mu, species, inc)

    # color_dict for stratification
    colors_dict = Dict("B1"=> "#BF42FC",
                       "B2"=> "#009E48",
                       "B3"=> "#E59201")
    c = [colors_dict[bName] for bName in bacNames]
    rC, gC, bC = HEX2RGB(c)
    colors_database = [RGB(rC[i]/255, gC[i]/255, bC[i]/255) for i in eachindex(rC)]
    colors = colors_database[species]
    alphas = [ac ? 1.0 : 0.2 for ac in active]

    plot()
    for i in eachindex(x)
        if !include_radius && !active[i]
            plot!(circleShape(x[i], y[i], 4.64e-7 * 1e6), seriestype =[:shape], aspect_rate=1, c=colors[i], linewidth=0.1, linealpha= 0.2, fillalpha=alphas[i], aspect_ratio=1, label=false)
        else
            plot!(circleShape(x[i], y[i], radius[i]), seriestype =[:shape], aspect_rate=1, c=colors[i], linewidth=0.1, fillalpha=alphas[i], aspect_ratio=1, label=false)
        end
    end

    for i in eachindex(bacNames)
        scatter!([],[], color = colors_database[i], linewidth=0.1, label=bacNames[i])
    end

    if !only_last
        plot!([],[], title="Time = $((i-1) * dT_save)", label=false)
    end

    fig = plot!(xlimits=xlim.*1e6, ylimits=ylim.*1e6, xlabel="Position along x-axis [µm]", ylabel="Position along y-axis [µm]",
                foreground_color_legend = nothing, background_color_legend = nothing, legendposition=:bottomright, dpi=600)

    directory = string(pwd(), @sprintf("\\results\\%04d", sim_number))
    filename = "$(directory)\\$(i).png"
    savefig(filename)
    return fig
end


function loaddata(sim_number, finished)

    directory = string(pwd(), @sprintf("\\results\\%04d\\results1D.jld2", sim_number))

    if finished
        sim_file = string(pwd(), @sprintf("\\results\\%04d\\sim_%04d.jld2", sim_number, sim_number))
    else
        sim_file = string(pwd(), @sprintf("\\sim_%04d.jld2", sim_number))
    end

    bac_vecint, bac_matfloat, bac_matint, bac_matbool = load(directory, "bac_saved_vecint", "bac_saved_matfloat", "bac_saved_matint", "bac_saved_matbool")
    grid_float, grid_int, constants_float, constants_vecstring = load(sim_file, "grid_float", "grid_int", "constants_float", "constants_vecstring")

    return bac_vecint, bac_matfloat, bac_matint, bac_matbool, grid_float, grid_int, constants_float, constants_vecstring

end


function getlimitdata(bac_vecint, bac_matfloat, grid_float, index)

    final_nBacs = bac_vecint.nBacs[index]
    println("Final number of bacteria: $(final_nBacs)")

    xlim = (minimum(bac_matfloat.x[index, 1:final_nBacs]) - 5*grid_float.dx,
            maximum(bac_matfloat.x[index, 1:final_nBacs]) + 5*grid_float.dx)
    ylim = (minimum(bac_matfloat.y[index, 1:final_nBacs]) - 5*grid_float.dy,
            maximum(bac_matfloat.y[index, 1:final_nBacs]) + 5*grid_float.dy)
    println(xlim, ylim)
    return xlim, ylim
end


function plotBacs(sim_number, finished, only_last, include_radius)
    """
    Function that plots the bacteria
    Arguments:
    sim_number              Integer indicating the simulation number
    finished                Boolean indicating whether simulation is finished (influences where data is taken from)
    only_last               Boolean whether only final state should be plotted, if false, a gif is made
    include_radius          Boolean indicating whether inactive bacteria should be plotted with their real radius (true) or a fixed radius (false). Fixed radius shows species of inactive cells better
    """

    bac_vecint, bac_matfloat, bac_matint, bac_matbool, grid_float, grid_int, constants_float, constants_vecstring = loaddata(sim_number, finished)
    
    lastnonzero = findlast(bac_vecint.nBacs .!= 0)
    xlim, ylim = getlimitdata(bac_vecint, bac_matfloat, grid_float, lastnonzero)

    if only_last
        fig = save_plot(lastnonzero, xlim, ylim, bac_vecint, bac_matfloat, bac_matint, bac_matbool, constants_vecstring.speciesNames, constants_float.dT_save, sim_number, only_last, include_radius)
    else
        # Gif implementation
        fig = "DONE"
        anim = @animate for i in 1:lastnonzero
            if bac_vecint.nBacs[i] != 0
                save_plot(i, xlim, ylim, bac_vecint, bac_matfloat, bac_matint, bac_matbool, constants_vecstring.speciesNames, constants_float.dT_save, sim_number, only_last, include_radius)
            end
        end
        gif(anim, string(pwd(), "\\results\\$(sim_number)\\$(sim_number).gif"), fps=10)

        for i in 1:lastnonzero
            rm(string(pwd(), "\\results\\$(sim_number)\\$(i).png"))
        end
    end

    println("DONE!")
    return fig
end
