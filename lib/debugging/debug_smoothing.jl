include(string(pwd(), "\\lib\\diffusion\\smoothing.jl"))
using BenchmarkTools

# L = [0 1 0; 1 -4 1; 0 1 0]
# # alpha = dT * diff_coef / (2 * dx^2)
# alpha = 2.65e-7 * 3.6e-6 / (2 * (2e-6)^2)
# L_0 = alpha .* L
# L_lhs = [0 0 0; 0 1 0; 0 0 0] .- ((1 ./ 2) .^(2 .* 0)) .* L_0

# phi = reshape(rand(131*131), 131,131)
# rhs = reshape(rand(131*131), 131,131)

# for i in 1:1000
#     smoothing(phi, rhs, L_lhs)
# end

using DSP
using ImageFiltering
a = [1 1 1 1; 1 2 2 1; 1 2 2 1; 1 1 1 1]
b = [0 1 0; 1 -4 1; 0 1 0]
@btime conv(a, b)[2:end-1, 2:end-1]
display(conv(a, b)[2:end-1, 2:end-1])
@btime @view conv(a, b)[2:end-1, 2:end-1]
display(@view conv(a, b)[2:end-1, 2:end-1])
@btime imfilter(a, reflect(centered(b)), Fill(0))
display(imfilter(a, reflect(centered(b)), Fill(0)))
@btime imfilter(a, reflect(centered(b)), Inner())
display(imfilter(a, reflect(centered(b)), Inner()))
