a = ones(4,4,2)
a[:,:,2] .= 2
b = [8, 6]
c = BitArray([0 0 0 0; 0 1 1 0; 1 0 0 1; 1 1 1 1])

special_b = reshape(b, 1,1,:)

special_c = repeat(c, outer = [1,1,2])

c .* a + .!special_c .* special_b
