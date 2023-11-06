using Test
import XLSX
loading_file = string(dirname(Base.source_dir()), "\\lib\\pre_processing\\loadPresetFile.jl")
include(loading_file)

create_mat_file = string(dirname(Base.source_dir()), "\\lib\\pre_processing\\create_mat.jl")
include(create_mat_file)

bacteria_shove_file = string(dirname(Base.source_dir()), "\\lib\\bacteria\\bacteria_shove.jl")
include(bacteria_shove_file)

struct General
    properties::Dict{Symbol, Any}
end
General() = General(Dict{Symbol, Any}())

Base.getproperty(x::General, property::Symbol) = getfield(x, :properties)[property]
Base.setproperty!(x::General, property::Symbol, value) = getfield(x, :properties)[property] = value
Base.propertynames(x::General) = keys(getfield(x, :properties))

filename = string(Base.source_dir(), "\\","test_file.xlsx")
grid, bac_init, constants, settings, init_params = loadPresetFile(filename)

@testset "rand_circle" begin
    n, x_centre, y_centre, radius = 30, 129, 129, 80
    candidate_x, candidate_y = rand_circle(n, x_centre, y_centre, radius) # Generate 30 points around centre 129,129 with max radius of 80
    @test size(candidate_x) == size(candidate_y) == (1, n) # Test size of Matrices
    @test typeof(candidate_x) == typeof(candidate_y) == Matrix{Float64} # Test type of returned Matrices
    @test all(sqrt.((candidate_x .- x_centre) .^2 + (candidate_y .- y_centre) .^2) .<= radius ) == true # Check all points within radius
end

@testset "blue_noise_circle" begin
    n, x_centre, y_centre, radius = 30, 129, 129, 80
    candidate_x, candidate_y = blue_noise_circle(n, x_centre, y_centre, radius) # Generate 30 points around centre 129,129 with max radius of 80
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
    new_bac = bacteria_shove(bac, grid, constants)
    @test size(new_bac.x) == size(new_bac.y)
    @test size(new_bac.x) == size(bac.x)
    @test typeof(new_bac.x) == typeof(new_bac.y)
    @test typeof(new_bac.x) == typeof(bac.x)
end

@testset "general" begin
    if settings.model_type == "suspension"
        @test size(bac.species) == (144,)
        @test size(bac.active) == (144,)
    else
        @test size(bac.species) == (118,)
        @test size(bac.active) == (118,)
    end
    @test all(bac.species .< 4) == true
    @test all(bac.active .== 1) == true
end
