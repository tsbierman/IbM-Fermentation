using Test
import XLSX
using InvertedIndices
using Random
using DSP

include(string(pwd(),"\\lib\\Lib_Module.jl"))
include(string(pwd(), "\\lib\\diffusion\\prolongation.jl"))

@testset "prolongation" begin
    phi_coarse = [1 5 3 6; 3 6 2 3; 4 2 2 3; 1 1 6 4]
    L_prolongation = [1 2 1; 2 4 2; 1 2 1] ./4
    sz = (7, 7)
    expected = [1 3 5 4 3 4.5 6; 2 3.75 5.5 4 2.5 3.5 4.5; 3 4.5 6 4 2 2.5 3; 3.5 3.75 4 3 2 2.5 3; 4 3 2 2 2 2.5 3; 2.5 2 1.5 2.75 4 3.75 3.5; 1 1 1 3.5 6 5 4]
    result = prolongation(phi_coarse, L_prolongation, sz)
    @test expected == round.(result, digits=2)
end
