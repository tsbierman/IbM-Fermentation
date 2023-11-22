using Test
import XLSX
using InvertedIndices
using Random

# Include the 4 detachment files

include(string(pwd(), "\\lib\\detachment\\computeRoot.jl"))
include(string(pwd(), "\\lib\\detachment\\recalculateT.jl"))
include(string(pwd(), "\\lib\\detachment\\calcTimeOfDetach.jl"))
include(string(pwd(), "\\lib\\detachment\\calculateLocalDetachmentRate.jl"))

create_mat_file = string(pwd(), "\\lib\\pre_processing\\create_mat.jl")
include(create_mat_file)

filename = string(Base.source_dir(), "\\","test_file.xlsx")
grid, bac, constants, settings, init_params = create_mat(filename)

@testset "computeRoot" begin
    @test(computeRoot(Inf, Inf, 1, 0.001) == Inf) # Artificial values for test
    @test(round(computeRoot(Inf, 10, 1, 0.001), digits=4) == 10.001)
    @test(round(computeRoot(10, Inf, 1, 0.001), digits=4) == 10.001)
    @test(round(computeRoot(10, 10, 1, 0.001), digits=4) == 10.0007)
end

@testset "calculateLocalDetachmentRate" begin
    @test calculateLocalDetachmentRate(30,30, constants.kDet, grid, 50.5*grid.dx, 50.5*grid.dx) != Inf # Not Inf
    @test calculateLocalDetachmentRate(30,51, constants.kDet, grid, 50.5*grid.dx, 50.5*grid.dx) > 0 # Always bigger than 0
    @test calculateLocalDetachmentRate(51,51, constants.kDet, grid, 50.5*grid.dx, 50.5*grid.dx) == Inf # Only centre will give Inf
    @test calculateLocalDetachmentRate(51.1,51, constants.kDet, grid, 50.5*grid.dx, 50.5*grid.dx) != Inf # Only centre will give Inf
end
