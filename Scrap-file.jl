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

