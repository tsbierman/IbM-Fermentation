using Test
import XLSX
using InvertedIndices
using Random
using DSP

include(string(pwd(), "\\inclusion_file.jl"))
include(string(pwd(),"\\lib\\Lib_Module.jl"))

create_mat_file = string(pwd(), "\\lib\\pre_processing\\create_mat.jl")
include(create_mat_file)

filename = string(pwd(), "\\test\\test_file.xlsx")
grid_float, grid_int, bac_vecfloat, bac_vecint, bac_vecbool, constants_float, constants_vecfloat, constants_vecint, constants_vecstring, constants_vecbool, constants_matfloat, settings_bool, settings_string, init_params = create_mat(filename, -1)
grid2bac, grid2nBacs = determine_where_bacteria_in_grid(grid_float, grid_int, bac_vecfloat)
diffusion_region, focus_region = determine_diffusion_region(grid2bac, grid2nBacs, bac_vecfloat, grid_float, grid_int)
Ks = constants_matfloat.Ks
Ki = constants_matfloat.Ki
Keq = constants_matfloat.Keq
chrM = constants_matfloat.chrM
chunks = 0
nChunks_dir = 4

Sh_ini = 1e-5
StV = [init_params.init_concs;1;0]
calculate_pH = settings_bool.pHincluded
Tol = constants_float.pHtolerance

@testset "determine_max_growth_rate_and_maint" begin
    res = determine_max_growth_rate_and_maint(1, 298, 1e-7) # Mock values for test
    @test length(res) == 2
    @test round(res[1], digits=4) == 0.0390
    @test round(res[2], digits=4) == 0.0082
end

@testset "calculate_monod" begin
    conc1 = [2e-3, 1e-3, 3e-3, 4e-3]
    conc1_1 = [2e-4, 1e-3, 3e-3, 4e-3]
    conc1_2 = [2e-5, 1e-3, 3e-3, 4e-3]
    conc3 = [2e-3, 0, 0, 0]
    Ks_test2 = deepcopy(Ks)
    Ki_test3 = deepcopy(Ki)
    Ks_test2[1,:] = [1e-4, 2e-5, 1e-4, 3e-4]
    Ki_test3[1,:] = [5e-3, 1e-2, 0, 0]
    @test round(calculate_monod(Ks[1,:], Ki[1,:], conc1), digits=3) == 0.995
    @test round(calculate_monod(Ks[1,:], Ki[1,:], conc1_1), digits=3) == 0.952
    @test round(calculate_monod(Ks[1,:], Ki[1,:], conc1_2), digits=3) == 0.667
    @test round(calculate_monod(Ks_test2[1,:], Ki[1,:], conc1_1), digits=3) == 0.588
    @test round(calculate_monod(Ks[1,:], Ki_test3[1,:], conc3), digits=3) == 0.711
    @test round(calculate_monod(Ks_test2[1,:], Ki_test3[1,:], conc1_1), digits=3) == 0.514
end

@testset "solve_pH" begin
    spcM = zeros(size(chrM))
    # StV = [0.0003333333333333333, 0.0003333333333333333, 0.0003333333333333333, 1.0e-20, 1.0e-20, 0.001, 0.0005, 0.001, 1.0, 0.0]
    StV = [1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1, 0]
    spcM_res = calculate_spcM!(spcM, 1e-7, Keq, StV)
    @test all(spcM_res .>=0)
    @test all(spcM_res[:,end] .== 0)
    @test all(spcM_res[1:5, 2] .== 1e-3)
    @test size(spcM_res) == (10, 5)
    @test all(spcM_res[1:5, Not(2)] .== 0)
    @test spcM_res[8,2] == 1e-3
    @test sum(spcM_res[6, :]) == 1e-3
    @test sum(spcM_res[7, :]) == 1e-3
    @test sum(spcM_res) == sum(StV)

    spcM, Sh = solve_pH(Sh_ini, StV, Keq, chrM, calculate_pH, Tol)
    @test 1 > Sh > 1e-14
end

@testset "calculate_reaction_matrix" begin
    pH = 3.0
    fake_init_conc = [1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3]
    fake_bulk_conc = [1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1, 0]
    fake_conc = ones(grid_int.ny, grid_int.nx, length(constants_vecstring.compoundNames)) * 1e-4
    fake_conc = set_concentrations!(fake_conc, fake_init_conc, diffusion_region)
    reaction_matrix, mu, pH = calculate_reaction_matrix!(grid2bac, grid2nBacs, bac_vecfloat, bac_vecint, bac_vecbool, diffusion_region, fake_conc, constants_float, constants_vecfloat, constants_vecint, constants_matfloat, pH, chunks, nChunks_dir, settings_bool)

    if settings_string.model_type in ("granule", "mature granule")
        @test size(reaction_matrix) == (257,257,8)
        @test length(mu) == length(bac_vecfloat.x)
        @test all(-1 .< mu .< 1)
        @test size(pH) == (257,257)
        @test all(3.0 .< pH .< 4.1)
        @test pH[1,1] > 3.45 # Should be the bulk 4.00
        @test pH[127,127] < 3.45 # Should be diffusion region 3.00
    end
end
