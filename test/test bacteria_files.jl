using Test
import XLSX
using InvertedIndices
using Random

include(string(pwd(),"\\lib\\bacteria\\Bacteria_Module.jl"))
include(string(pwd(), "\\lib\\pre_processing\\create_mat.jl"))

filename = string(pwd(), "\\test\\test_file.xlsx")

@testset "bacteria_divide" begin
    grid_float, grid_int, bac_vecfloat, bac_vecint, bac_vecbool, constants_float, constants_vecfloat, constants_vecint, constants_vecstring, constants_vecbool, constants_matfloat, settings_bool, settings_string, init_params = create_mat(filename, -1)

    bac_vecfloat.molarMass[1] = constants_float.max_bac_mass_grams / constants_float.bac_MW * 1.1
    bac_vecfloat.radius[1] = ( (bac_vecfloat.molarMass[1] * constants_float.bac_MW / constants_float.bac_rho) * (3 / (4 * pi))) ^ (1/3)
    bac_vecfloat.mu = ones(length(bac_vecfloat.x)) * 0.02

    old_bac_vecfloat = deepcopy(bac_vecfloat)
    old_bac_vecint = deepcopy(bac_vecint)
    old_bac_vecbool = deepcopy(bac_vecbool)

    bac_vecfloat, bac_vecint, bac_vecbool, cycle = Bacteria_Module.bacteria_divide!(bac_vecfloat, bac_vecint, bac_vecbool, constants_float)

    @test length(bac_vecfloat.x) == length(old_bac_vecfloat.x) + 1
    @test cycle == 1
    @test bac_vecint.species[1] == old_bac_vecint.species[1] == bac_vecint.species[end]
    @test 0.45 < bac_vecfloat.molarMass[1] / old_bac_vecfloat.molarMass[1] < 0.55
    @test 0.45 < bac_vecfloat.molarMass[end] / old_bac_vecfloat.molarMass[1] < 0.55
    @test (bac_vecfloat.radius[1] < old_bac_vecfloat.radius[1]) & (bac_vecfloat.radius[end] < old_bac_vecfloat.radius[1])

    grid_float, grid_int, bac_vecfloat, bac_vecint, bac_vecbool, constants_float, constants_vecfloat, constants_vecint, constants_vecstring, constants_vecbool, constants_matfloat, settings_bool, settings_string, init_params = create_mat(filename, -1)
    bac_vecfloat.molarMass[1] = constants_float.max_bac_mass_grams / constants_float.bac_MW * 1.1
    bac_vecfloat.molarMass[2] = constants_float.max_bac_mass_grams / constants_float.bac_MW * 1.1
    bac_vecfloat.molarMass[3] = constants_float.max_bac_mass_grams / constants_float.bac_MW * 1.1
    bac_vecfloat.mu = ones(length(bac_vecfloat.x)) * 0.02
    old_bac_vecfloat = deepcopy(bac_vecfloat)
    bac_vecfloat, bac_vecint, bac_vecbool, cycle = Bacteria_Module.bacteria_divide!(bac_vecfloat, bac_vecint, bac_vecbool, constants_float)
    @test length(bac_vecfloat.x) == length(old_bac_vecfloat.x) + 3

end

@testset "killBacs" begin
    grid_float, grid_int, bac_vecfloat, bac_vecint, bac_vecbool, constants_float, constants_vecfloat, constants_vecint, constants_vecstring, constants_vecbool, constants_matfloat, settings_bool, settings_string, init_params = create_mat(filename, -1)

    old_bac_vecfloat = deepcopy(bac_vecfloat)
    my_zeros = BitArray(zeros(length(bac_vecfloat.x)))
    my_zeros[1] = 1
    bac_vecfloat.mu = ones(length(bac_vecfloat.x)) * 0.02
    bac_vecfloat, bac_vecint, bac_vecbool = Bacteria_Module.killBacs!(bac_vecfloat, bac_vecint, bac_vecbool, my_zeros)
    @test length(bac_vecfloat.x) == length(old_bac_vecfloat.x) - 1
    @test bac_vecfloat.x[1] == old_bac_vecfloat.x[2]

    grid_float, grid_int, bac_vecfloat, bac_vecint, bac_vecbool, constants_float, constants_vecfloat, constants_vecint, constants_vecstring, constants_vecbool, constants_matfloat, settings_bool, settings_string, init_params = create_mat(filename, -1)

    old_bac_vecfloat = deepcopy(bac_vecfloat)
    indices = [1;5;10]
    bac_vecfloat.mu = ones(length(bac_vecfloat.x)) * 0.02
    bac_vecfloat, bac_vecint, bac_vecbool = Bacteria_Module.killBacs!(bac_vecfloat, bac_vecint, bac_vecbool, indices)
    @test length(bac_vecfloat.x) == length(old_bac_vecfloat.x) - 3
    @test bac_vecfloat.x[1] == old_bac_vecfloat.x[2]
end

@testset "bacteria die" begin
    grid_float, grid_int, bac_vecfloat, bac_vecint, bac_vecbool, constants_float, constants_vecfloat, constants_vecint, constants_vecstring, constants_vecbool, constants_matfloat, settings_bool, settings_string, init_params = create_mat(filename, -1)

    bac_vecfloat.molarMass[1] = constants_float.min_bac_mass_grams / constants_float.bac_MW * 0.9
    old_bac_vecfloat = deepcopy(bac_vecfloat)
    bac_vecfloat.mu = ones(length(bac_vecfloat.x)) * 0.02
    bac_vecfloat, bac_vecint, bac_vecbool = Bacteria_Module.bacteria_die!(bac_vecfloat, bac_vecint, bac_vecbool, constants_float)
    @test length(old_bac_vecfloat.x) - 1 == length(bac_vecfloat.x)
    @test old_bac_vecfloat.x[2] == bac_vecfloat.x[1]

    grid_float, grid_int, bac_vecfloat, bac_vecint, bac_vecbool, constants_float, constants_vecfloat, constants_vecint, constants_vecstring, constants_vecbool, constants_matfloat, settings_bool, settings_string, init_params = create_mat(filename, -1)

    bac_vecfloat.molarMass[1] = constants_float.min_bac_mass_grams / constants_float.bac_MW * 0.9
    bac_vecfloat.molarMass[10] = constants_float.min_bac_mass_grams / constants_float.bac_MW * 0.9
    bac_vecfloat.molarMass[20] = constants_float.min_bac_mass_grams / constants_float.bac_MW * 0.9
    old_bac_vecfloat = deepcopy(bac_vecfloat)
    bac_vecfloat.mu = ones(length(bac_vecfloat.x)) * 0.02
    bac_vecfloat, bac_vecint, bac_vecbool = Bacteria_Module.bacteria_die!(bac_vecfloat, bac_vecint, bac_vecbool, constants_float)
    @test length(old_bac_vecfloat.x) - 3 == length(bac_vecfloat.x)
end

@testset "bacteria_inactivate" begin
    grid_float, grid_int, bac_vecfloat, bac_vecint, bac_vecbool, constants_float, constants_vecfloat, constants_vecint, constants_vecstring, constants_vecbool, constants_matfloat, settings_bool, settings_string, init_params = create_mat(filename, -1)

    bac_vecfloat.mu = ones(length(bac_vecfloat.x)) * 0.02
    bac_vecfloat.molarMass[1] = constants_float.min_bac_mass_grams / constants_float.bac_MW * 0.9
    bac_vecbool.active[2] = false
    bac_vecfloat.mu[2] = -0.01
    bac_vecbool = Bacteria_Module.bacteria_inactivate!(bac_vecfloat, bac_vecbool, constants_float)
    @test bac_vecbool.active[1] == 0
    @test bac_vecbool.active[2] == 0
end

@testset "update_bacterial" begin
    grid_float, grid_int, bac_vecfloat, bac_vecint, bac_vecbool, constants_float, constants_vecfloat, constants_vecint, constants_vecstring, constants_vecbool, constants_matfloat, settings_bool, settings_string, init_params = create_mat(filename, -1)

    bac_vecfloat.mu = ones(length(bac_vecfloat.x)) * 0.02
    bac_vecfloat.mu[2] = -0.02
    old_bac_vecfloat = deepcopy(bac_vecfloat)
    bac_vecfloat = Bacteria_Module.update_bacterial_mass!(bac_vecfloat, bac_vecbool, 0.05)
    bac_vecfloat = Bacteria_Module.update_bacterial_radius!(bac_vecfloat, constants_float)
    @test bac_vecfloat.molarMass[1] > old_bac_vecfloat.molarMass[1]
    @test bac_vecfloat.molarMass[2] < old_bac_vecfloat.molarMass[2]
    @test bac_vecfloat.radius[1] > old_bac_vecfloat.radius[1]
    @test bac_vecfloat.radius[2] < old_bac_vecfloat.radius[2]
end

@testset "bacteria_detachment" begin
    grid_float, grid_int, bac_vecfloat, bac_vecint, bac_vecbool, constants_float, constants_vecfloat, constants_vecint, constants_vecstring, constants_vecbool, constants_matfloat, settings_bool, settings_string, init_params = create_mat(filename, -1)
    bac_vecfloat.mu = ones(length(bac_vecfloat.x)) * 0.2
    old_bac_vecfloat = deepcopy(bac_vecfloat)
    if settings_string.detachment in ("SBR", "none")
        bac_vecfloat, bac_vecint, bac_vecbool = Bacteria_Module.bacteria_detachment!(bac_vecfloat, bac_vecint, bac_vecbool, grid_float, grid_int, constants_float, settings_string, constants_float.dT_bac, init_params.invHRT[1])
        @test bac_vecfloat.x == old_bac_vecfloat.x
    end

    if settings_string.detachment == "naive"
        bac_vecfloat.x[1] = 1e-6
        bac_vecfloat, bac_vecint, bac_vecbool = Bacteria_Module.bacteria_detachment!(bac_vecfloat, bac_vecint, bac_vecbool, grid_float, grid_int, constants_float, settings_string, constants_float.dT_bac, init_params.invHRT[1])
        @test length(bac_vecfloat.x) == length(old_bac_vecfloat.x) - 1
        @test bac_vecfloat.x[1] == old_bac_vecfloat.x[2]
    end

    if settings_string.detachment == "suspension"
        bac_vecfloat.mu[1] = 0.02
        bac_vecfloat, bac_vecint, bac_vecbool = Bacteria_Module.bacteria_detachment!(bac_vecfloat, bac_vecint, bac_vecbool, grid_float, grid_int, constants_float, settings_string, constants_float.dT_bac, init_params.invHRT[1])
        @test length(bac_vecfloat.x) == length(old_bac_vecfloat.x) - 1
        @test bac_vecfloat.x[1] == old_bac_vecfloat.x[2]
    end
end
