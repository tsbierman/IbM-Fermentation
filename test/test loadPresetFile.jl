using Test
import XLSX
loading_file = string(dirname(Base.source_dir()), "\\lib\\pre-processing\\loadPresetFile.jl")
include(loading_file)

struct General
    properties::Dict{Symbol, Any}
end
General() = General(Dict{Symbol, Any}())

Base.getproperty(x::General, property::Symbol) = getfield(x, :properties)[property]
Base.setproperty!(x::General, property::Symbol, value) = getfield(x, :properties)[property] = value
Base.propertynames(x::General) = keys(getfield(x, :properties))


@testset "loadPresetFile" begin
    filename = string(Base.source_dir(), "\\","test_file.xlsx")
    grid, bac_init, constants, settings, init_params = loadPresetFile(filename)

    #Test Discretization
    @test grid.dx == 2e-6
    @test round(constants.max_granule_radius[1], digits=6)  == 2.53e-4
    @test round(constants.dT, digits=9) == 2.65e-7
    @test round(constants.dynamicDT.maxDT, digits=9) == 2.12e-7

    # Test diffusion
    @test round(constants.diffusion_rates[5], digits=9) == 7.56e-6
    @test length(constants.compoundNames) == 8

    # Test Parameters
    @test constants.T == 293.15
    @test init_params.invHRT == 0.1

    # Test Bacteria
    @test constants.kDist == 1

    # Test solver
    @test isnan(constants.pHtolerance) == true
    @test settings.type == "Neut"

    # Test influent and initial condition
    @test constants.influent_concentrations[7] == 5e-4
    @test init_params.init_bulk_conc[7] == 5e-4

    # Test thermodynamic parameters
    @test round(constants.Keq[6,1]) == 5.88e2
    @test constants.Keq[7,2] == 100
    @test constants.chrM[1,1] == 0
    @test constants.chrM[6,4] == -2
    
    # Test Ks
    @test constants.Ks[1,1] == 1.00e-5
    @test constants.Ks[2,2] == 1.00e-5

    # Test matrix
    @test constants.MatrixMet[1,1] == -100
    @test size(constants.MatrixMet) == (10,3)
end
