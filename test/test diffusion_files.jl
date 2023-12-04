using Test
import XLSX
using InvertedIndices
using Random
using DSP

include(string(pwd(),"\\lib\\Lib_Module.jl"))
include(string(pwd(), "\\lib\\diffusion\\prolongation.jl"))
include(string(pwd(), "\\lib\\diffusion\\restriction.jl"))

@testset "prolongation" begin
    phi_coarse = [1 5 3 6; 3 6 2 3; 4 2 2 3; 1 1 6 4]
    L_prolongation = [1 2 1; 2 4 2; 1 2 1] ./4
    sz = (7, 7)
    expected = [1 3 5 4 3 4.5 6; 2 3.75 5.5 4 2.5 3.5 4.5; 3 4.5 6 4 2 2.5 3; 3.5 3.75 4 3 2 2.5 3; 4 3 2 2 2 2.5 3; 2.5 2 1.5 2.75 4 3.75 3.5; 1 1 1 3.5 6 5 4]
    result = prolongation(phi_coarse, L_prolongation, sz)
    @test expected == round.(result, digits=2)
end

@testset "restriction" begin
    phi_fine = [1 2 5 1 1; 3 2 3 4 2; 4 4 4 5 1; 5 3 3 2 5; 1 5 2 4 3]
    L_restriction = [1 2 1; 2 4 2; 1 2 1] ./ 16
    expected = [16 38 14; 45 57 34; 27 37 32] ./ 16
    result = restriction(phi_fine, L_restriction)
    @test expected == round.(result, digits=4)
end
