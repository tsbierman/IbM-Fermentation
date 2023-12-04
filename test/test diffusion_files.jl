using Test
import XLSX
using InvertedIndices
using Random
using DSP

include(string(pwd(),"\\lib\\Lib_Module.jl"))
include(string(pwd(), "\\lib\\diffusion\\prolongation.jl"))
include(string(pwd(), "\\lib\\diffusion\\restriction.jl"))
include(string(pwd(), "\\lib\\diffusion\\residual.jl"))
include(string(pwd(), "\\lib\\diffusion\\smoothing.jl"))
include(string(pwd(), "\\lib\\diffusion\\create_dirichlet_boundary.jl"))

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

@testset "create_dirichlet_boundary" begin
    mock_matrix = ones(3,3)
    mock_matrix1 = zeros(50,31)
    value = 8
    result = create_dirichlet_boundary(mock_matrix, value)
    result1 = create_dirichlet_boundary(mock_matrix1, value)
    @test size(result) == (5,5)
    @test sum(result) == sum(mock_matrix) + (2 * size(mock_matrix, 1) + 2 * size(mock_matrix, 2) + 4) * value
    @test all(result[2:end-1, 2:end-1] .== 1)
    @test all(result[1,:] .== value)
    @test all(result[:,1] .== value)
    @test all(result[end,:] .== value)
    @test all(result[:,end] .== value)

    @test size(result1) == (52,33)
    @test sum(result1) == sum(mock_matrix1) + (2 * size(mock_matrix1, 1) + 2 * size(mock_matrix1, 2) + 4) * value
    @test all(result1[2:end-1, 2:end-1] .== 0)
    @test all(result1[1,:] .== value)
    @test all(result1[:,1] .== value)
    @test all(result1[end,:] .== value)
    @test all(result1[:,end] .== value)
end

@testset "residual" begin
    phi = [1 5 3 6; 3 6 2 3; 4 2 2 3; 1 1 6 4] # Mock concentrations
    L_lhs = [1 2 1; 2 4 2; 1 2 1] / 4          # Mock stencil
    rhs = [3 5 2 9; 8 4 2 1; 0 3 2 5; 6 5 4 1]
    expected_result = [-3.5 -6.25 -9.75 -0.5; -2.25 -10.5 -11 -8.75; -8.75 -8.5 -10 -4.5; 2 -2 -6.75 -8]
    result = residual(phi, rhs, L_lhs)
    @test round.(result,digits=2) == expected_result
end
