# file_loc = string(Base.source_dir(), "\\","test_file.xlsx")
# grid, bac_init, constants, settings, init_params = loadPresetFile(file_loc);

# a = reshape([11;2;13;14;15], 5, :)
# b = reshape([10;20;12;4;20], :, 5)

# d = zeros(4,1)
# println(typeof(d[1:end]))
# c = d[1:end] .- b

# d = findmax(minimum(c, dims = 1))
# a = transpose(range(2,11,10))
# b = transpose(range(11,20,10))

using InvertedIndices
a = range(1,15,15)
# b = reshape(a, 3,5)

c = a .>= 7
a = a[Not(c)]
# a

# rand((1:5), (15,1))

# file_loc = string(Base.source_dir(), "\\planning\\test_file.xlsx")
# grid, bac_init, constants, settings, init_params = loadPresetFile(file_loc)

# findall(constants.speciesNames .== "B2")[1]

# typeof(length(constants.speciesNames))
# using Random
# a = randperm(15)
# b = [1;15;3;8;4;5]
# print(size(b))
# a = reshape(a, 15,1)
# println(a)
# I = sortperm(b)
# a[I]
# zeros(size(a))
# x, y = rand_circle(20, 30, 30, 15)
# d[1] = 30
# d[1:end] .- x
# println(size(x))

# X = zeros(30, 1)
# X[1] = 20

# print(typeof(X))

# a = range(1, 15, 15)
# a = reshape(a, 3, 5)
# sqrt.(a)


# a = randperm(15)
# a = reshape(a, 15,:)
# dx = rand(1:30)
# dy = rand(1:30)
# println(dx,"\n", dy)
# option_a = (dx * dy)
# option_b = ((dx/2)^2 + (dy/2)^2)/2 * pi
# option_c = ((dx + dy)/2 / 2)^2 * pi
# println(option_a,"\n", option_b,"\n", option_c)



# a = General()

# a.number = 15

# function changer(a)
#     a.number = 16
# end

# changer(a)
# a.number
