root_folder = string(pwd(), "\\testing_push_walk_map")
subdirectories = Set(Base.Iterators.flatten(collect(walkdir(root_folder))))

for subdir in subdirectories
    push!(LOAD_PATH, string(subdir))
end

a = 2
b = 3
c = 5

println(addition_func(a, b, c))
println(multi_func(a, b))
