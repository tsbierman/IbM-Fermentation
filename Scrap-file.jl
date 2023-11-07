a = ones(15)
a[5] = 10
a[9] = 2
b = a .> 1
d = a .> 5
c = b .& d

typeof(b)
e = BitArray(ones(15))
