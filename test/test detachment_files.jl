using Test
import XLSX
using InvertedIndices
using Random

# Include the 4 detachment files

include(string(pwd(), "\\lib\\detachment\\computeRoot.jl"))
include(string(pwd(), "\\lib\\detachment\\recalculateT.jl"))
include(string(pwd(), "\\lib\\detachment\\calcTimeOfDetach.jl"))
include(string(pwd(), "\\lib\\detachment\\calculateLocalDetachmentRate.jl"))

include(string(pwd(), "\\lib\\pre_processing\\create_mat.jl"))

filename = string(pwd(), "\\test\\test_file.xlsx")
grid_float, grid_int, bac_vecfloat, bac_vecint, bac_vecbool, constants_float, constants_vecfloat, constants_vecint, constants_vecstring, constants_vecbool, constants_matfloat, settings_bool, settings_string, init_params = create_mat(filename, -1)

@testset "computeRoot" begin
    @test(computeRoot(Inf, Inf, 1, 0.001) == Inf) # Artificial values for test
    @test(round(computeRoot(Inf, 10, 1, 0.001), digits=4) == 10.001)
    @test(round(computeRoot(10, Inf, 1, 0.001), digits=4) == 10.001)
    @test(round(computeRoot(10, 10, 1, 0.001), digits=4) == 10.0007)
end

@testset "calculateLocalDetachmentRate" begin
    @test calculateLocalDetachmentRate(30,30, constants_float.kDet, grid_float, 50.5*grid_float.dx, 50.5*grid_float.dx) != Inf # Not Inf
    @test calculateLocalDetachmentRate(30,51, constants_float.kDet, grid_float, 50.5*grid_float.dx, 50.5*grid_float.dx) > 0 # Always bigger than 0
    @test calculateLocalDetachmentRate(51,51, constants_float.kDet, grid_float, 50.5*grid_float.dx, 50.5*grid_float.dx) == Inf # Only centre will give Inf
    @test calculateLocalDetachmentRate(51.1,51, constants_float.kDet, grid_float, 50.5*grid_float.dx, 50.5*grid_float.dx) != Inf # Only centre will give Inf
end
