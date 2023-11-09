# a = rand(5,3) * 2
# display(a)

# b = findall((a .> 0.5) .& (a .< 1.5))
# b = 0.5 .< a .< Inf
# for i in eachindex(b)
#     println(i)
#     println("DONE")
# end

# b[3]

# a = zeros(3,3,3)
# a[1,1,2] = 1
# a[1,1,3] = 10
# c = findall(a[1,1,:] .!= 0)
# b = a[1,1,c]

# for value in b
#     println(value)
# end

2 + 1:2 + 3

a = 1:15
c = [1;3;4;5;8;10;14;15]
b = BitArray([1;0;1;1;1;0;0;1;0;1;0;0;0;1;1])
using InvertedIndices

a[Not(c)]
mean(a)


a = BitArray([1;1;1;0;0])
b = BitArray([0;0;1;1;1])

a .& .!b
