using Test
import XLSX
using InvertedIndices
using Random
using DSP

include(string(pwd(),"\\lib\\Lib_Module.jl"))

create_mat_file = string(dirname(Base.source_dir()), "\\lib\\pre_processing\\create_mat.jl")
include(create_mat_file)

filename = string(Base.source_dir(), "\\","test_file.xlsx")
grid, bac, constants, settings, init_params = create_mat(filename)
grid2bac, grid2nBacs = Lib_Module.determine_where_bacteria_in_grid(grid, bac)
diffusion_region, focus_region = Lib_Module.determine_diffusion_region(grid2bac, grid2nBacs, bac, grid)

@testset "set_concentrations" begin
    conc = zeros(grid.ny, grid.nx, length(constants.compoundNames))
    conc1 = Lib_Module.set_concentrations!(conc, init_params.init_concs, diffusion_region)
    conc2 = Lib_Module.set_concentrations!(conc, init_params.init_bulk_conc, .!diffusion_region)
    conc3 = Lib_Module.set_concentrations!(conc, init_params.init_concs, BitArray(zeros(grid.ny, grid.nx)))

    @test round(conc1[Int(ceil(grid.ny/2)), Int(ceil(grid.nx/2)), 1],digits=6) == 3.33e-4
    @test round(conc1[Int(ceil(grid.ny/2)), Int(ceil(grid.nx/2)), 6],digits=6) == 1.00e-3
    @test all(conc1[1,1,:] .== 0)

    @test round(conc2[1,1, 1],digits=6) == 3.33e-4
    @test round(conc2[1,1, 6],digits=6) == 1.00e-3
    @test all(conc2[Int(ceil(grid.ny/2)), Int(ceil(grid.nx/2)), :] .== 0)
    @test conc3 == conc
end

@testset "calculate_bulk_concentrations" begin # Can only test for the initial case as reactionMatrix = 0 in that case (calculate_reactionMatrix not yet implemented)
    bulk_concs, invHRT = Lib_Module.calculate_bulk_concentrations(bac, constants, init_params.init_bulk_conc, init_params.invHRT,0, constants.dT_bac, settings)
    @test bulk_concs == init_params.init_bulk_conc
    @test all(bulk_concs .> 0)
end

@testset "determine_where_bacteria_in_grid" begin
    total_bacs = sum(grid2nBacs)
    total_non_zeros = sum(grid2bac .!= 0)
    @test size(grid2bac) == (grid.ny, grid.nx, 4) # Not 4 if this is adapted
    @test size(grid2nBacs) == (grid.ny, grid.nx)
    @test total_bacs == total_non_zeros == length(bac.x)
end

@testset "determine_diffusion_region" begin
    @test size(diffusion_region) == (grid.ny, grid.nx)
    @test sum(diffusion_region) > length(bac.x)/4 # max 4 per grid cell
    
    @test focus_region.x0 < (ceil(minimum(bac.x) / grid.dx))
    @test focus_region.x1 > (ceil(maximum(bac.x) / grid.dx))
    @test focus_region.y0 < (ceil(minimum(bac.y) / grid.dy))
    @test focus_region.y1 > (ceil(maximum(bac.y) / grid.dy))
end

@testset "create_chunks" begin
    chunks = Lib_Module.create_chunks(10, focus_region)
    print(typeof(chunks))
    @test typeof(chunks) == Main.Lib_Module.General
    @test size(chunks.indices_x) == size(chunks.indices_x) == (10,2)
    @test length(chunks.dx_chunk) == length(chunks.dy_chunk) == 1
end
