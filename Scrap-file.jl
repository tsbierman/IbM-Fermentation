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
