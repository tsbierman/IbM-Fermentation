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

# a = zeros(10, 2)

# for i in axes(a)[1]
#     println(i)
# end

# b = 1:4
# c = reshape(b, 1,1,:)
# d = BitArray(ones(10,10))
# d[5,5] = 0
# e = .!d .* c

# a = 1:10
# println(size(a))
# b = [1:10]
# println(size(b))
# for i in b
#     println(i)
# end

# a = range(1,50,50)
# b = reshape(a, 10, 5)
# c = b.> 8

# kernel = ones(3,3) ./ -8
# kernel[2,2] = 1
# # println(b)
# # b[2:end-1, 2:end-1]

# using DSP
# conv(c, kernel)[2:end-1, 2:end-1] .> 1e-15

# d = BitArray(ones(5,3))
# d[1:6] .= 0
# findall(d .!= 0)

# a = [1, 2, 3, 4, 5, 6]

# b = repeat(a, 6)
# c = transpose(repeat(a, 6))


# x = transpose(repeat(a, 1, 6))
# y = transpose(repeat(a, inner=(1,6)))

# tempX = nodeMidCoordinatesX[1:5]
# x = repeat(tempX,inner=4)
# x = repeat(tempX,4)
# println(x)

# a = BitArray(ones(3,3))
# a[[1 3 4 7]] .= 0
# println(a)
# b = transpose(a)
# collect(Base.Iterators.flatten(b))

# a = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18]
# b = reshape(a, 2, 3, 3)
# c = sum(b, dims=(1,2))
# dropdims(c, dims=(1,2))

a = ones(4,5) * 2
b = ones(4,5)

b[[3 4 5 10]] .= 5
c = a .* b
sum(c)