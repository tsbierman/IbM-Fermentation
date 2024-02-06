using Test
import XLSX
using InvertedIndices
using Random

include(string(pwd(),"\\lib\\bacteria\\Bacteria_Module.jl"))

include(string(pwd(), "\\lib\\pre_processing\\create_mat.jl"))

filename = string(pwd(), "\\test\\test_file.xlsx")
@testset "bacteria_divide" begin
    grid, bac, constants, settings, init_params = create_mat(filename, -1)
    bac.molarMass[1] = constants.max_bac_mass_grams/constants.bac_MW * 1.1
    bac.radius[1] = ( (bac.molarMass[1] * constants.bac_MW / constants.bac_rho) * (3 / (4 * pi))) ^ (1/3)
    bac.mu = ones(length(bac.x)) * 0.02
    old_bac = deepcopy(bac)
    bac, cycle = Bacteria_Module.bacteria_divide!(bac, constants)
    @test length(bac.x) == length(old_bac.x) + 1
    @test cycle == 1
    @test bac.species[1] == old_bac.species[1] == bac.species[end]
    @test 0.45 < bac.molarMass[1]/old_bac.molarMass[1] < 0.55
    @test 0.45 < bac.molarMass[end]/old_bac.molarMass[1] < 0.55
    @test (bac.radius[1] < old_bac.radius[1]) & (bac.radius[end] < old_bac.radius[1])

    grid, bac, constants, settings, init_params = create_mat(filename, -1)
    bac.molarMass[1] = constants.max_bac_mass_grams/constants.bac_MW * 1.1
    bac.molarMass[2] = constants.max_bac_mass_grams/constants.bac_MW * 1.1
    bac.molarMass[3] = constants.max_bac_mass_grams/constants.bac_MW * 1.1
    bac.mu = ones(length(bac.x)) * 0.02
    old_bac = deepcopy(bac)
    bac, cycle = Bacteria_Module.bacteria_divide!(bac, constants)
    @test length(bac.x) == length(old_bac.x) + 3

end

@testset "killBacs" begin
    grid, bac, constants, settings, init_params = create_mat(filename, -1)
    old_bac = deepcopy(bac)
    my_zeros = BitArray(zeros(length(bac.x)))
    my_zeros[1] = 1
    bac.mu = ones(length(bac.x)) * 0.02
    bac = Bacteria_Module.killBacs!(bac, my_zeros)
    @test length(bac.x) == length(old_bac.x) - 1
    @test bac.x[1] == old_bac.x[2]

    grid, bac, constants, settings, init_params = create_mat(filename, -1)
    old_bac = deepcopy(bac)
    indices = [1;5;10]
    bac.mu = ones(length(bac.x)) * 0.02
    bac = Bacteria_Module.killBacs!(bac, indices)
    @test length(bac.x) == length(old_bac.x) - 3
    @test bac.x[1] == old_bac.x[2]
end

@testset "bacteria die" begin
    grid, bac, constants, settings, init_params = create_mat(filename, -1)
    bac.molarMass[1] = constants.min_bac_mass_grams / constants.bac_MW * 0.9
    old_bac = deepcopy(bac)
    bac.mu = ones(length(bac.x)) * 0.02
    bac = Bacteria_Module.bacteria_die!(bac, constants)
    @test length(old_bac.x) - 1 == length(bac.x)
    @test old_bac.x[2] == bac.x[1]

    grid, bac, constants, settings, init_params = create_mat(filename, -1)
    bac.molarMass[1] = constants.min_bac_mass_grams / constants.bac_MW * 0.9
    bac.molarMass[10] = constants.min_bac_mass_grams / constants.bac_MW * 0.9
    bac.molarMass[20] = constants.min_bac_mass_grams / constants.bac_MW * 0.9
    old_bac = deepcopy(bac)
    bac.mu = ones(length(bac.x)) * 0.02
    bac = Bacteria_Module.bacteria_die!(bac, constants)
    @test length(old_bac.x) - 3 == length(bac.x)
end

@testset "bacteria_inactivate" begin
    grid, bac, constants, settings, init_params = create_mat(filename, -1)
    bac.mu = ones(length(bac.x)) * 0.02
    bac.molarMass[1] = constants.min_bac_mass_grams / constants.bac_MW * 0.9
    bac.active[2] = false
    bac.mu[2] = -0.01
    bac = Bacteria_Module.bacteria_inactivate!(bac, constants)
    @test bac.active[1] == 0
    @test bac.active[2] == 0
end

@testset "update_bacterial" begin
    grid, bac, constants, settings, init_params = create_mat(filename, -1)
    bac.mu = ones(length(bac.x)) * 0.02
    bac.mu[2] = -0.02
    old_bac = deepcopy(bac)
    bac = Bacteria_Module.update_bacterial_mass!(bac, 0.05)
    bac = Bacteria_Module.update_bacterial_radius!(bac, constants)
    @test bac.molarMass[1] > old_bac.molarMass[1]
    @test bac.molarMass[2] < old_bac.molarMass[2]
    @test bac.radius[1] > old_bac.radius[1]
    @test bac.radius[2] < old_bac.radius[2]
end

@testset "bacteria_detachment" begin
    grid, bac, constants, settings, init_params = create_mat(filename, -1)
    bac.mu = ones(length(bac.x)) * 0.2
    old_bac = deepcopy(bac)
    if settings.detachment in ("SBR", "none")
        bac = Bacteria_Module.bacteria_detachment!(bac, grid, constants, settings, constants.dT_bac, init_params.invHRT)
        @test bac.x == old_bac.x
    end

    if settings.detachment == "naive"
        bac.x[1] = 1e-6
        bac = Bacteria_Module.bacteria_detachment!(bac, grid, constants, settings, constants.dT_bac, init_params.invHRT)
        @test length(bac.x) == length(old_bac.x) - 1
        @test bac.x[1] == old_bac.x[2]
    end

    if settings.detachment == "suspension"
        bac.mu[1] = 0.02
        bac = Bacteria_Module.bacteria_detachment!(bac, grid, constants, settings, constants.dT_bac, init_params.invHRT)
        @test length(bac.x) == length(old_bac.x) - 1
        @test bac.x[1] == old_bac.x[2]
    end
end
