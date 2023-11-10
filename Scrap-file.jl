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

a = zeros(10, 2)

for i in axes(a)[1]
    println(i)
end
