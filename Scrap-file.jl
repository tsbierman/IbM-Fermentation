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
