using Plots
using StatsPlots

# This document plots the data from literature to visualise ti
BO_growth_rates = [0.069, 0.032, 0.07518, 0.0358, 0.032]
AM_growth_rates = [0.05, 0.058, 0.028, 0.058, 0.01475, 0.0141667]
HM_growth_rates = [0.126, 0.03, 0.05, 0.068, 0.051, 0.069, 0.1675, 0.17]

# Names and colors for Anaerobic Digestion
species = ["BO", "AM", "HM"]
color_list = ["#A020F0" "#E8BB00"  "#3bd6c6"]

# Calculate the chosen numbers
BO_chosen = 1.2/24
AM_chosen = 0.4/24
HM_chosen = 2.1/24

# Horizontal shift to prevent overlap
BO_shift = (rand(size(BO_growth_rates, 1)) .* 2 .- 1) .* 0.05
AM_shift = (rand(size(AM_growth_rates, 1)) .* 2 .- 1) .* 0.05
HM_shift = (rand(size(HM_growth_rates, 1)) .* 2 .- 1) .* 0.05

# Plot the points in a scatter plot, with some horizontal noise
scatter(ones(size(BO_growth_rates, 1)) .+ BO_shift, BO_growth_rates, color = color_list[1], label=nothing)
scatter!(ones(size(AM_growth_rates, 1)).*2 .+ AM_shift, AM_growth_rates, color = color_list[2], label=nothing, xlim=(0.5,3.5), ylim=(0,0.2))
scatter!(ones(size(HM_growth_rates, 1)).*3 .+ HM_shift, HM_growth_rates, color = color_list[3], label=nothing, xlabel="Specie", ylabel="Maximum growth rate (h-1)", xticks=(1:3, species))
plot!([0.90, 1.1], [BO_chosen, BO_chosen], linestyle=:dash, color = :black, linewidth=2, label=nothing)
plot!([1.90, 2.1], [AM_chosen, AM_chosen], linestyle=:dash, color = :black, linewidth=2, label=nothing)
plot!([2.90, 3.1], [HM_chosen, HM_chosen], linestyle=:dash, color = :black, linewidth=2, label=nothing)

savefig(string(pwd(), "\\growth_rate_literature.png"))
