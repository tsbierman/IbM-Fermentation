using Test
import XLSX
using InvertedIndices
using Random
using DSP

include(string(pwd(),"\\lib\\Lib_Module.jl"))
include(string(pwd(), "\\inclusion_file.jl"))
filename = string(pwd(), "\\test\\test_file.xlsx")
grid_float, grid_int, bac_vecfloat, bac_vecint, bac_vecbool, constants_float, constants_vecfloat, constants_vecint, constants_vecstring, constants_vecbool, constants_matfloat, settings_bool, settings_string, init_params = create_mat(filename, -1)

# debug_struct = General()
# constants.debug = debug_struct
# constants.debug.plotDiffRegion = false

grid2bac, grid2nBacs = Lib_Module.determine_where_bacteria_in_grid(grid_float, grid_int, bac_vecfloat)
diffusion_region, focus_region = Lib_Module.determine_diffusion_region(grid2bac, grid2nBacs, bac_vecfloat, grid_float, grid_int)

@testset "set_concentrations" begin
    if settings_string.model_type in ("granule", "mature granule")
        concs = zeros(grid_int.ny, grid_int.nx, length(constants_vecstring.compoundNames))
        conc1 = Lib_Module.set_concentrations!(concs, init_params.init_concs, diffusion_region)
        conc2 = Lib_Module.set_concentrations!(concs, init_params.init_bulk_conc, .!diffusion_region)
        conc3 = Lib_Module.set_concentrations!(concs, init_params.init_concs, BitArray(zeros(grid_int.ny, grid_int.nx)))

        @test round(conc1[Int(ceil(grid_int.ny/2)), Int(ceil(grid_int.nx/2)), 1],digits=6) == 3.33e-4
        @test round(conc1[Int(ceil(grid_int.ny/2)), Int(ceil(grid_int.nx/2)), 6],digits=6) == 1.00e-3
        @test all(conc1[1,1,:] .== 0)

        @test round(conc2[1,1, 1],digits=6) == 3.33e-4
        @test round(conc2[1,1, 6],digits=6) == 1.00e-3
        @test all(conc2[Int(ceil(grid_int.ny/2)), Int(ceil(grid_int.nx/2)), :] .== 0)
        @test conc3 == concs
    end
end

@testset "calculate_bulk_concentrations" begin # Can only test for the initial case as reactionMatrix = 0 in that case (calculate_reactionMatrix not yet implemented)
    bulk_conc, invHRT = Lib_Module.calculate_bulk_concentrations(bac_vecfloat, bac_vecbool, constants_float, constants_vecfloat, constants_vecint, constants_vecstring, constants_vecbool, constants_matfloat, init_params.init_bulk_conc, init_params.invHRT[1], 0, constants_float.dT_bac, settings_bool, settings_string)
    @test bulk_conc == init_params.init_bulk_conc
    @test all(bulk_conc .> 0)
end

@testset "determine_where_bacteria_in_grid" begin
    total_bacs = sum(grid2nBacs)
    total_non_zeros = sum(grid2bac .!= 0)
    @test size(grid2bac) == (grid_int.ny, grid_int.nx, 9) # Not 9 if this is adapted
    @test size(grid2nBacs) == (grid_int.ny, grid_int.nx)
    @test total_bacs == total_non_zeros == length(bac_vecfloat.x)
end

@testset "determine_diffusion_region" begin
    @test size(diffusion_region) == (grid_int.ny, grid_int.nx)
    @test sum(diffusion_region) > length(bac_vecfloat.x)/4 # max 4 per grid cell
    
    @test focus_region.x0 < (ceil(minimum(bac_vecfloat.x) / grid_float.dx))
    @test focus_region.x1 > (ceil(maximum(bac_vecfloat.x) / grid_float.dx))
    @test focus_region.y0 < (ceil(minimum(bac_vecfloat.y) / grid_float.dy))
    @test focus_region.y1 > (ceil(maximum(bac_vecfloat.y) / grid_float.dy))
end

@testset "create_chunks" begin
    chunks_matrix, chunks_int = Lib_Module.create_chunks(10, focus_region)
    @test size(chunks_matrix.indices_x) == size(chunks_matrix.indices_x) == (10,2)
    @test length(chunks_int.dx_chunk) == length(chunks_int.dy_chunk) == 1
end

conc, bulk_concs, invHRT, reaction_matrix, pH, bac_vecfloat, bac_vecint, bac_vecbool = initTime!(grid_float, grid_int, bac_vecfloat, bac_vecint, bac_vecbool, init_params, constants_float, constants_vecfloat, constants_vecint, constants_vecstring, constants_vecbool, constants_matfloat, settings_bool, settings_string)

@testset "initTime" begin
    @test size(conc) == (grid_int.ny, grid_int.nx, length(constants_vecstring.compoundNames))
    @test conc[1,1,1] == bulk_concs[1]
    @test round(conc[257,257,1],digits=6) == 3.33e-4

    @test length(bulk_concs) == length(constants_vecstring.compoundNames)
    @test bulk_concs[1] == init_params.init_bulk_conc[1]

    @test length(invHRT) == 1
    
    @test size(reaction_matrix) == size(conc)

    @test size(pH) == (grid_int.ny, grid_int.nx)
    @test pH[1,1] == constants_float.pHsetpoint

    @test length(bac_vecfloat.mu) == length(bac_vecfloat.x)
end
