using Test
import XLSX
using InvertedIndices
using Random
using DSP
using Printf

include(string(pwd(), "\\inclusion_file.jl"))
Temp_Time = Float_struct()

filename = string(pwd(), "\\test\\test_file.xlsx")
grid_float, grid_int, bac_vecfloat, bac_vecint, bac_vecbool, constants_float, constants_vecfloat, constants_vecint, constants_vecstring, constants_vecbool, constants_matfloat, settings_bool, settings_string, init_params = create_mat(filename, -1)

@testset "bulk_conc_diff_within_limit" begin
    conc1 = 100 
    conc2 = 103
    conc3 = 101
    conc4 = 98.1
    temp_vector = zeros(3)
    temp_vector[:] = [conc2, conc3, conc4]
    temp_vector2 = zeros(2)
    temp_vector2[:] = [conc3, conc4]

    duo_standard = zeros(2)
    duo_standard[:] = [conc1, conc1]
    trio_standard = zeros(3)
    trio_standard[:] = [conc1, conc1, conc1]

    @test bulk_conc_diff_within_limit(conc2, conc1, constants_float) == false
    @test bulk_conc_diff_within_limit(conc3, conc1, constants_float) == true
    @test bulk_conc_diff_within_limit(conc4, conc1, constants_float) == true
    @test bulk_conc_diff_within_limit(temp_vector, trio_standard, constants_float) == false
    @test bulk_conc_diff_within_limit(temp_vector2, duo_standard, constants_float) == true
end

@testset "decrease_dT_bac" begin
    Temp_Time.current = 15
    Temp_Time.changed_dT_bac = 10
    Temp_Time.minDT_bac = 0.1
    Temp_Time.dT_bac = 0.5

    decrease_dT_bac!(Temp_Time, "Testing enough space for time change")
    @test Temp_Time.dT_bac == 0.4
    @test Temp_Time.changed_dT_bac == 15

    Temp_Time.current = 16
    Temp_Time.dT_bac = 0.12
    decrease_dT_bac!(Temp_Time, "Testing not enough space for time change")
    @test Temp_Time.dT_bac == 0.1
    @test Temp_Time.changed_dT_bac == 16
end

@testset "decrease_dT_diffusion" begin
    Temp_Time.current = 15
    Temp_Time.changed_dT = 10
    Temp_Time.minDT = 5.29e-9
    Temp_Time.dT = 1e-8

    decrease_dT_diffusion!(Temp_Time, "Testing enough space for time change", grid_float.dx, constants_vecfloat)
    @test round(Temp_Time.dT,digits=9) == 0.9e-8
    @test Temp_Time.changed_dT == 15

    Temp_Time.current = 16
    Temp_Time.dT = 5.4e-9
    decrease_dT_diffusion!(Temp_Time, "Testing not enough space for time change", grid_float.dx, constants_vecfloat)
    @test Temp_Time.dT == 5.29e-9
    @test Temp_Time.changed_dT == 16

    Temp_Time.current = 17
    Temp_Time.dT = 5.29e-9
    decrease_dT_diffusion!(Temp_Time, "Testing already at minimum", grid_float.dx, constants_vecfloat)
    @test Temp_Time.dT == 5.29e-9
    @test Temp_Time.changed_dT == 16
end

@testset "increase_dT_bac" begin
    Temp_Time.current = 15
    Temp_Time.changed_dT_bac = 10
    Temp_Time.maxDT_bac = 1
    Temp_Time.dT_bac = 0.5

    increase_dT_bac!(Temp_Time, "Testing enough space for time change")
    @test Temp_Time.dT_bac == 0.625
    @test Temp_Time.changed_dT_bac == 15

    Temp_Time.current = 16
    Temp_Time.dT_bac = 0.81
    increase_dT_bac!(Temp_Time, "Testing not enough space for time change")
    @test Temp_Time.dT_bac == 1
    @test Temp_Time.changed_dT_bac == 16
end

@testset "increase_dT_diffusion" begin
    Temp_Time.current = 15
    Temp_Time.changed_dT = 10
    Temp_Time.maxDT = 2.12e-7
    Temp_Time.dT = 1e-8

    increase_dT_diffusion!(Temp_Time, "Testing enough space for time change", grid_float.dx, constants_vecfloat)
    @test round(Temp_Time.dT,digits=10) == 1.11e-8
    @test Temp_Time.changed_dT == 15

    Temp_Time.current = 16
    Temp_Time.dT = 2e-7
    increase_dT_diffusion!(Temp_Time, "Testing not enough space for time change", grid_float.dx, constants_vecfloat)
    @test Temp_Time.dT == 2.12e-7
    @test Temp_Time.changed_dT == 16
end

@testset "multiple_high_iters" begin
    Temp_Time.dT_bac = 0.5
    Temp_Time.current = 15
    Temp_Time.changed_dT = 10
    iProf1 = 5
    iProf2 = 2
    iDiffusion1 = 25
    iDiffusion2 = 15
    nDiffIters1 = [5, 8, 23, 25,]
    nDiffIters2 = [30, 50, 19, 25,]
    nDiffIters3 = [30, 30, 30, 30,]

    @test multiple_high_iters(iDiffusion1, iProf1, nDiffIters1, Temp_Time, constants_vecint) == true
    @test multiple_high_iters(iDiffusion2, iProf1, nDiffIters3, Temp_Time, constants_vecint) == false
    @test multiple_high_iters(iDiffusion1, iProf1, nDiffIters2, Temp_Time, constants_vecint) == false
    @test multiple_high_iters(iDiffusion1, iProf2, nDiffIters1, Temp_Time, constants_vecint) == false
    Temp_Time.changed_dT = 14
    @test multiple_high_iters(iDiffusion1, iProf1, nDiffIters1, Temp_Time, constants_vecint) == false
end

@testset "multiple_low_initRES" begin
    Temp_Time.dT_bac = 0.5
    Temp_Time.current = 15
    Temp_Time.changed_dT_bac = 10
    iProf1 = 3
    iProf2 = 2
    maxInitRES1 = [0.1, 0.1, 0.1,]
    maxInitRES2 = [0.1, 0.2, 0.1,]
    maxInitRES3 = [0.1, 0.1, 0.1,]

    @test multiple_low_initRES(iProf1, maxInitRES1, Temp_Time, constants_float, constants_vecint) == true
    @test multiple_low_initRES(iProf1, maxInitRES2, Temp_Time, constants_float, constants_vecint) == false
    @test multiple_low_initRES(iProf2, maxInitRES1, Temp_Time, constants_float, constants_vecint) == false
    Temp_Time.changed_dT_bac = 14
    @test multiple_low_initRES(iProf1, maxInitRES1, Temp_Time, constants_float, constants_vecint) == false
end

@testset "non_convergent" begin
    tol = 1e-7
    iRES5 = 5
    iRES6 = 6
    RES_both_fail = ones(3,6)
    RES_both_fail[2,6] = 10
    RES_both_fail[1,5] = 10 + 1e-5
    RES_both_fail[3,4] = 10 + 2e-7
    @test non_convergent(iRES6, RES_both_fail, tol) == false

    RES_pov1_succ = deepcopy(RES_both_fail)
    RES_pov1_succ[1,5] = 10 + 0.9e-7
    @test non_convergent(iRES6, RES_pov1_succ, tol) == false

    RES_pov2_succ = deepcopy(RES_both_fail)
    RES_pov2_succ[3,4] = 10 + 1.9e-7
    @test non_convergent(iRES6, RES_pov2_succ, tol) == false

    RES_pov_succeed = deepcopy(RES_both_fail)
    RES_pov_succeed[1,5] = 10 + 0.9e-7
    RES_pov_succeed[3,4] = 10 + 1.9e-7
    @test non_convergent(iRES6, RES_pov_succeed, tol) == true

    RES_neg1_succ = deepcopy(RES_both_fail)
    RES_neg1_succ[1,5] = 10 - 1e-6
    @test non_convergent(iRES6, RES_neg1_succ, tol) == false

    RES_neg2_succ = deepcopy(RES_both_fail)
    RES_neg1_succ[1,5] = 10 + 1e-6
    @test non_convergent(iRES6, RES_neg2_succ, tol) == false

    RES_neg_succeed = deepcopy(RES_both_fail)
    RES_neg_succeed[1,5] = 10 - 1e-4
    @test non_convergent(iRES6, RES_neg_succeed, tol) == true

    RES_succeed = deepcopy(RES_both_fail)
    RES_succeed[1,5] = 10 - 0.9e-7
    RES_succeed[3,4] = 10 + 1.9e-7
    @test non_convergent(iRES6, RES_succeed, tol) == true
    @test non_convergent(iRES5, RES_succeed, tol) == false
end

@testset "slow_convergence" begin
    RES_allsucc = ones(3,250)
    RES_allsucc[1,250] = 10
    RES_allsucc[3,249] = 10 + 1e-4
    RES_allsucc[3,248] = 10 + 2e-4
    iRES250 = 250
    iRES200 = 200
    @test slow_convergence(iRES250, RES_allsucc, constants_float, constants_vecint) == true
    @test slow_convergence(iRES200, RES_allsucc, constants_float, constants_vecint) == false
    RES_little_fail = deepcopy(RES_allsucc)
    RES_little_fail[3,248] = 10 + 2e-2
    @test slow_convergence(iRES250, RES_little_fail, constants_float, constants_vecint) == false
    RES_direc_fail = deepcopy(RES_allsucc)
    RES_direc_fail[3,248] = 10 - 2e-2
    @test slow_convergence(iRES250, RES_direc_fail, constants_float, constants_vecint) == false
    RES_no_conv_true = deepcopy(RES_allsucc)
    RES_no_conv_true[3,249] = 10 + 0.9e-7
    RES_no_conv_true[3,248] = 10 + 1.9e-7
    @test slow_convergence(iRES250, RES_no_conv_true, constants_float, constants_vecint) == false
end

@testset "upward_trend" begin
    iRES3 = 3
    iRES2 = 2
    RES = ones(3,3)
    @test upward_trend(iRES3, RES) == false
    RES[2,3] = 1.2
    @test upward_trend(iRES2, RES) == false
    @test upward_trend(iRES3, RES) == true
end
