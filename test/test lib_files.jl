using Test
import XLSX
using InvertedIndices
using Random
using DSP

include(string(pwd(),"\\lib\\Lib_Module.jl"))

include(string(pwd(), "\\lib\\initTime.jl"))
include(string(pwd(), "\\lib\\reaction_matrix\\calculate_reaction_matrix.jl"))
include(string(pwd(), "\\lib\\reaction_matrix\\solve_pH.jl"))
include(string(pwd(), "\\lib\\reaction_matrix\\rMatrix_section.jl"))
include(string(pwd(), "\\lib\\reaction_matrix\\calculate_monod.jl"))
include(string(pwd(), "\\lib\\reaction_matrix\\determine_max_growth_rate_and_maint.jl"))

include(string(pwd(), "\\lib\\pre_processing\\create_mat.jl"))

filename = string(pwd(), "\\test\\test_file.xlsx")
grid, bac, constants, settings, init_params = create_mat(filename)

debug_struct = General()
constants.debug = debug_struct
constants.debug.plotDiffRegion = false

grid2bac, grid2nBacs = Lib_Module.determine_where_bacteria_in_grid(grid, bac)
diffusion_region, focus_region = Lib_Module.determine_diffusion_region(grid2bac, grid2nBacs, bac, grid)

@testset "set_concentrations" begin
    concs = zeros(grid.ny, grid.nx, length(constants.compoundNames))
    conc1 = Lib_Module.set_concentrations!(concs, init_params.init_concs, diffusion_region)
    conc2 = Lib_Module.set_concentrations!(concs, init_params.init_bulk_conc, .!diffusion_region)
    conc3 = Lib_Module.set_concentrations!(concs, init_params.init_concs, BitArray(zeros(grid.ny, grid.nx)))

    @test round(conc1[Int(ceil(grid.ny/2)), Int(ceil(grid.nx/2)), 1],digits=6) == 3.33e-4
    @test round(conc1[Int(ceil(grid.ny/2)), Int(ceil(grid.nx/2)), 6],digits=6) == 1.00e-3
    @test all(conc1[1,1,:] .== 0)

    @test round(conc2[1,1, 1],digits=6) == 3.33e-4
    @test round(conc2[1,1, 6],digits=6) == 1.00e-3
    @test all(conc2[Int(ceil(grid.ny/2)), Int(ceil(grid.nx/2)), :] .== 0)
    @test conc3 == concs
end

@testset "calculate_bulk_concentrations" begin # Can only test for the initial case as reactionMatrix = 0 in that case (calculate_reactionMatrix not yet implemented)
    bulk_conc, invHRT = Lib_Module.calculate_bulk_concentrations(bac, constants, init_params.init_bulk_conc, init_params.invHRT,0, constants.dT_bac, settings)
    @test bulk_conc == init_params.init_bulk_conc
    @test all(bulk_conc .> 0)
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

conc, bulk_concs, invHRT, reaction_matrix, pH, bac = initTime!(grid, bac, init_params, constants, settings)

@testset "initTime" begin
    @test size(conc) == (grid.ny, grid.nx, length(constants.compoundNames))
    @test conc[1,1,1] == bulk_concs[1]
    @test round(conc[257,257,1],digits=6) == 3.33e-4

    @test length(bulk_concs) == length(constants.compoundNames)
    @test bulk_concs[1] == init_params.init_bulk_conc[1]

    @test length(invHRT) == 1
    
    @test size(reaction_matrix) == size(conc)

    @test size(pH) == (grid.ny, grid.nx)
    @test pH[1,1] == constants.pHsetpoint

    @test length(bac.mu) == length(bac.x)
end
