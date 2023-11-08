b = randperm(15)
println(b)

c = sortperm(b)
println(c)
b[c]
