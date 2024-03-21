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
