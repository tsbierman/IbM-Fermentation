# ---------------------------- START PLOT HELP -----------------------------------
# a = [1, 2, 3, 4, 5, 6]

# b = repeat(a, 6)
# c = transpose(repeat(a, 6))


# x = transpose(repeat(a, 1, 6))
# y = transpose(repeat(a, inner=(1,6)))

# tempX = nodeMidCoordinatesX[1:5]
# x = repeat(tempX,inner=4)
# x = repeat(tempX,4)
# println(x)
# ---------------------------- END PLOT HELP -----------------------------------

# -------------------------- START PRIORITYQUEUE --------------------------
# using DataStructures
# pq = PriorityQueue()

# for i in range(10,1,10)
#     push!(pq, string(Int(i+1), ",", Int(i+1)) => i)
# end
# push!(pq, "4,89" => 4.90)
# pq[string(4, ",", 89)] = 50
# pq

# offSet = [0  0  1 -1; 1 -1 0 0]
# for nb in 1:length(offSet[1,:])
#     x_index_nb = 5 + offSet[1, nb]
#     y_index_nb = 5 + offSet[2, nb]

#     println("Index x $(x_index_nb)")
#     println("Index y $(y_index_nb)")

# end

# c = findfirst(a, ",")

# d = "236,5"
# e = findfirst(",", d)
# f = parse(Int64, d[1:e[1]-1])
# g = parse(Int64, d[e[1]+1:end])
# println(typeof(f))
# --------------------------END PRIORITYQUEUE----------------------------

# --------------START MULTIDIMENSIONAL BITARRAYS (FOR AFTER CONVOLUTION) -----------------
# a = ones(4,4,2)
# a[:,:,2] .= 2
# b = [8, 6]
# c = BitArray([0 0 0 0; 0 1 1 0; 1 0 0 1; 1 1 1 1])

# special_b = reshape(b, 1,1,:)

# special_c = repeat(c, outer = [1,1,2])

# c .* a + .!special_c .* special_b

# --------------END MULTIDIMENSIONAL BITARRAYS (FOR AFTER CONVOLUTION) -----------------

# ------------------------------- START SAVING AND LOADING -----------------------------

# using FileIO
# include(string(pwd(), "\\lib\\Struct_Module.jl"))

# a = 15
# b = [18, 31, 72]
# c = General()
# c.first = 1
# c.last = 18
# c.helppls = ones(4,4)
# c.helppls[2:3,2:3] .= 10

# save("testing_save.jld2", "a", a, "b", b, "c", c)
# # save("testing_save.jld2", "d", 1000) Extra save overwrites original
# a, b, c = load("testing_save.jld2", "a", "b", "c")

# println("DONE")

# ------------------------------- START SAVING AND LOADING -----------------------------
