using Test
import XLSX

include(string(pwd(), "\\lib\\pre_processing\\loadPresetFile.jl"))
include(string(pwd(), "\\lib\\Struct_Module.jl"))

@testset "loadPresetFile" begin
    filename = string(pwd(), "\\test\\test_file.xlsx")
    grid_float, grid_int, bac_init_float, bac_init_int, constants_float, constants_vecfloat, constants_vecint, constants_vecstring, constants_vecbool, constants_matfloat, settings_bool, settings_string, init_params = loadPresetFile(filename)

    #Test Discretization
    @test grid_float.dx == 2e-6
    @test round(constants_float.max_granule_radius, digits=6)  == 2.53e-4
    @test round(constants_float.dT, digits=9) == 2.65e-7
    @test round(constants_float.maxDT, digits=9) == 2.12e-7

    # Test diffusion
    @test round(constants_vecfloat.diffusion_rates[5], digits=9) == 7.56e-6
    @test length(constants_vecstring.compoundNames) == 8

    # Test Parameters
    @test constants_float.T == 293.15
    @test init_params.invHRT[1] == 0.1

    # Test Bacteria
    @test constants_float.kDist == 1

    # Test solver
    if settings_bool.pHincluded
        @test isnan(constants_float.pHtolerance) == false
    else
        @test isnan(constants_float.pHtolerance) == true
    end
    if settings_bool.structure_model
        @test settings_string.type == "Neut"
    end
        
    # Test influent and initial condition
    @test constants_vecfloat.influent_concentrations[7] == 5e-4
    @test init_params.init_bulk_conc[7] == 5e-4

    # Test thermodynamic parameters
    @test round(constants_matfloat.Keq[6,1]) == 5.88e2
    @test constants_matfloat.Keq[7,2] == 100
    @test constants_matfloat.chrM[1,1] == 0
    @test constants_matfloat.chrM[6,4] == -2
    
    # Test Ks
    @test constants_matfloat.Ks[1,1] == 1.00e-5
    @test constants_matfloat.Ks[2,2] == 1.00e-5

    # Test matrix
    @test constants_matfloat.MatrixMet[1,1] == -100
    @test size(constants_matfloat.MatrixMet) == (8,3)
end
