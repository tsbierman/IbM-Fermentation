using Test
import XLSX

include(string(pwd(), "\\inclusion_file.jl"))
filename = string(pwd(), "\\test\\test_file.xlsx")
grid_float, grid_int, bac_vecfloat, bac_vecint, bac_vecbool, constants_float, constants_vecfloat, constants_vecint, constants_vecstring, constants_vecbool, constants_matfloat, settings_bool, settings_string, init_params = create_mat(filename, -1)

@testset "rand_circle" begin
    n, x_centre, y_centre, radius = 30, 129, 129, 80
    candidate_x, candidate_y = rand_circle(n, x_centre, y_centre, radius) # Generate 30 points around centre 129,129 with max radius of 80
    @test size(candidate_x) == size(candidate_y) == (1, n) # Test size of Matrices
    @test typeof(candidate_x) == typeof(candidate_y) == Matrix{Float64} # Test type of returned Matrices
    @test all(sqrt.((candidate_x .- x_centre) .^2 + (candidate_y .- y_centre) .^2) .<= radius ) == true # Check all points within radius
end

@testset "blue_noise_circle" begin
    n, x_centre, y_centre, radius = 30, 129, 129, 80
    candidate_x, candidate_y = blue_noise_circle(n, x_centre, y_centre, radius) # Generate 30 points around centre (129,129) with max radius of 80
    @test size(candidate_x) == size(candidate_y) == (n,) # Test size of Vectors
    @test typeof(candidate_x) == typeof(candidate_y) == Vector{Float64} # Test type of returned Vectors
    @test all(sqrt.((candidate_x .- x_centre) .^2 + (candidate_y .- y_centre) .^2) .<= radius ) == true # Check all points within radius
    @test candidate_x[1] == x_centre
end


@testset "distribute_microcolonies" begin
    nColonies, nBacPerCol, r_colony, xrange, yrange = 8, 13, 1.2, [30, 220], [30, 220]
    x, y = distribute_microcolonies(nColonies, nBacPerCol, r_colony, xrange, yrange)
    @test size(x) == size(y) == (nColonies * nBacPerCol, ) # Test size of Vectors
    @test typeof(x) == typeof(y) == Vector{Float64} # Test type of returned Vectors
    @test all(xrange[1]-30 .<= x .<= xrange[2]+30) == true # Cannot test exactly due to noise, thus check whether on grid.
    @test all(yrange[1]-30 .<= y .<= yrange[2]+30) == true

    nColonies, nBacPerCol, r_colony, xrange, yrange = 1, 13, 1.2, [30, 220], [30,220]
    x, y = distribute_microcolonies(nColonies, nBacPerCol, r_colony, xrange, yrange)
    @test all(sqrt.((x .- x[1]) .^2 + (y .- y[1]) .^2) .<= r_colony ) == true # Check all points within radius
end

@testset "bacteria_shove" begin
    new_bac_vecfloat = bacteria_shove!(bac_vecfloat, grid_float, grid_int, constants_float)
    @test size(new_bac_vecfloat.x) == size(new_bac_vecfloat.y)
    @test size(new_bac_vecfloat.x) == size(bac_vecfloat.x)
    @test typeof(new_bac_vecfloat.x) == typeof(new_bac_vecfloat.y)
    @test typeof(new_bac_vecfloat.x) == typeof(bac_vecfloat.x)
end

@testset "general" begin
    if settings_string.model_type == "suspension"
        @test length(bac_vecint.species) <= 144
        @test length(bac_vecbool.active) <= 144
    else
        @test length(bac_vecint.species) <= 118
        @test length(bac_vecbool.active) <= 118
    end
    @test all(bac_vecint.species .< 4) == true
    @test all(bac_vecbool.active .== 1) == true
end
