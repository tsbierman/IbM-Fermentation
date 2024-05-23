using BenchmarkTools
# a = ["string1", "string2", "string3"]
# c = ["string3", "string4", "string5"]

# @time c[a .== "string2"][1]
# @time c[findfirst(==("string2"), a)]

# a = [1, 2, 3,4,5,6]
# b = [2,2,2,2,2,2]
# @time a * 5 ./ b
# @time a .* 5 ./ b


# function func1(x)
#     A = zeros(x,x)
#     for i in 1:x
#         for j in 1:x
#             A[i,j] = 1
#         end
#     end
#     return A
# end

# function func2(x)
#     A = zeros(x,x)
#     for i in 1:x
#         for j in 1:x
#             A[j, i] = 1
#         end
#     end
#     return A
# end

# @btime func1(10000);
# @btime func2(10000);

# a = [1 0 -1; 1 0 0; -1 0.5 0.5]
# b = [1,2,3]

# c = a .* b
# d = sum(c, dims=1)
# e = dropdims(sum(c, dims=1), dims=1)
# display(c)
# display(d)
# display(e)

# f = [1,2]
# g = zeros(size(b))
# g[1:length(f)] = f
# display(g)

# file = string(pwd(), "\\results\\9965\\results1D.jld2")
# bac_saved_vecint = load(file, "bac_saved_vecint")
# bac_saved_vecint.nBacs

# a = [1,2,3,4,5,10,0,0]
# minimum(a)
# findlast(a .!= 0)
# a .!= 0

# using Plots
# plot()
# x = [1,3,2]
# y = [1,3,2]

# function HEX2RGB(c)

#     rC, gC, bC = zeros(length(c)), zeros(length(c)), zeros(length(c))
#     for idx in eachindex(c)
#         color = c[idx]
#         RGB = parse(UInt, color[2:end], base=16)
#         rC[idx] = (RGB >> 16) & 0xFF
#         gC[idx] = (RGB >> 8) & 0xFF
#         bC[idx] = RGB & 0xFF
#     end
#     return rC, gC, bC
# end

# a = ["#BF42FC","#009E48","#E59201"]
# rC, gC, bC = HEX2RGB(a)
# colors = [RGB(rC[i]/255, gC[i]/255, bC[i]/255) for i in eachindex(rC)]

# color1 = RGB(1,0,0)
# color2 = RGB(0,1,0)
# color3 = RGB(0,0,1)
# color_match = [color1, color2, color3]
# indices = [1,2,3]
# println(color_match[indices])
# scatter!(x, y, c=colors[indices])
# plot!()

# a = [1,0,1,0,0,1,1,0]
# b = [i == 1 ? 5 : 4 for i in a]


# function circleShape(h, k, r)
#     """
#     Function to plot circles
#     """
#     thet = LinRange(0, 2*pi, 500)
#     return h .+ r * sin.(thet), k .+ r*cos.(thet)
# end

# using Plots
# plot()
# x_coor = [1,2,3]
# y_coor = [1,2,3]
# r = [1,1,1]
# for i in eachindex(x_coor)
#     plot!(circleShape(x_coor[i], y_coor[i], r[i]), seriestype = [:shape,], aspect_ratio = 1)
# end
# plot!()

# a = [1,2,3,4,5,6,7,8,9]
# bool1 = BitArray([1,1,1,1,1,0,0,0,0])
# bool2 = [-1,1,0,1,0,0,-1,1,1]
# a[bool1 .& (bool2 .!= 1)]

# using Plots

# R = 8.314/1000
# T = 298.15
# but_base = 48.3
# pro_base = 76.1
# meth_base = -135.926

# log_range = 1:7
# num_range = [10.0^(-num) for num in log_range]

# but = [but_base + R * T * log(num^2) for num in num_range]
# pro = [pro_base + R * T * log(num^3) for num in num_range]
# meth = [meth_base + R * T * log(num^(-4)) for num in num_range]

# hline([0], linestyle=:dash, c=:black, label=false)
# plot!(log_range, but, label="Butyrate oxidation", linewidth=3)
# plot!(log_range, pro, label="Propionate oxidation", linewidth=3)
# plot!(log_range, meth, label="Methanogenesis", linewidth=3)
# xflip!(true)
# plot!(xlabel="-log10(pH2)", ylabel="kJ/reaction", legend_position=:bottomleft)
# savefig("Gibbs energy of formation.png")

# variable = zeros(3,3)
# display(variable)
# variable[4,4] = 5
# display(variable)

# variable1 = zeros(3,3)
# variable2 = copy(variable1)
# display(variable2)
# variable1[2,2] = 5
# display(variable1)
# display(variable2)
# variable2[3,3] = 4
# display(variable1)
# display(variable2)

# x = [1,2,3,4,5,6]
# y1 = 2 * x .+ 15
# y2 = 2 .* x .+ 15

# species = [1,2,3,2,3,1,2,3,1,2]
# active = BitArray([1,0,0,1,1,1,1,1,0,1])
# radius = [5,8,9,7,5,6,3,2,5,8]
# active_species = active .& (species .== 1)
# mass_per_species = [sum(radius[active .& (species .== i)] * 2) for i in 1:maximum(species)]
# active_per_species = [sum(active .& (species .== i)) for i in 1:maximum(species)]

# a = ["a", "b", "c", "d"]
# c = [1 2 3 4 5; 6 7 8 9 10]

# b = findall(a .== "b")[1]
# c[b,3]

# function my_func1()
#     x = [1,5,7,8,9,6,5]
#     centre = [3]
#     pos = [1,2,3,2,3,2,1]
#     dist = zeros(length(x))
#     for index in 1:maximum(pos)
#         loca = pos .== index
#         dist[loca] = x[loca] .- centre[index]
#     end
#     return dist
# end

# function my_func2()
#     x = [1,5,7,8,9,6,5]
#     centre = [3.0]
#     dist = x .- centre
#     return dist
# end

# @btime(my_func1)
# @btime(my_func2)
