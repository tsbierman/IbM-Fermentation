using FileIO
include(string(pwd(), "\\lib\\Struct_Module.jl"))

a = 15
b = [18, 31, 72]
c = General()
c.first = 1
c.last = 18
c.helppls = ones(4,4)
c.helppls[2:3,2:3] .= 10

save("testing_save.jld2", "a", a, "b", b, "c", c)
# save("testing_save.jld2", "d", 1000) Extra save overwrites original
a, b, c = load("testing_save.jld2", "a", "b", "c")

println("DONE")
