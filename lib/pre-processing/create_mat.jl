#Initialise structs that will have to be used later
# Structs need to be declared at top level, so as this
struct Grid
    properties::Dict{Symbol, Any}
end
Grid() = Grid(Dict{Symbol, Any}())

Base.getproperty(x::Grid, property::Symbol) = getfield(x, :properties)[property]
Base.setproperty!(x::Grid, property::Symbol, value) = getfield(x, :properties)[property] = value
Base.propertynames(x::Grid) = keys(getfield(x, :properties))

struct Constants
    properties::Dict{Symbol, Any}
end
Constants() = Constants(Dict{Symbol, Any}())

Base.getproperty(x::Constants, property::Symbol) = getfield(x, :properties)[property]
Base.setproperty!(x::Constants, property::Symbol, value) = getfield(x, :properties)[property] = value
Base.propertynames(x::Constants) = keys(getfield(x, :properties))

struct Settings
    properties::Dict{Symbol, Any}
end
Settings() = Settings(Dict{Symbol, Any}())

Base.getproperty(x::Settings, property::Symbol) = getfield(x, :properties)[property]
Base.setproperty!(x::Settings, property::Symbol, value) = getfield(x, :properties)[property] = value
Base.propertynames(x::Settings) = keys(getfield(x, :properties))

struct Init_params
    properties::Dict{Symbol, Any}
end
Init_params() = Init_params(Dict{Symbol, Any}())

Base.getproperty(x::Init_params, property::Symbol) = getfield(x, :properties)[property]
Base.setproperty!(x::Init_params, property::Symbol, value) = getfield(x, :properties)[property] = value
Base.propertynames(x::Init_params) = keys(getfield(x, :properties))

struct Bac_init
    properties::Dict{Symbol, Any}
end
Bac_init() = Bac_init(Dict{Symbol, Any}())

Base.getproperty(x::Bac_init, property::Symbol) = getfield(x, :properties)[property]
Base.setproperty!(x::Bac_init, property::Symbol, value) = getfield(x, :properties)[property] = value
Base.propertynames(x::Bac_init) = keys(getfield(x, :properties))

struct Bac
    properties::Dict{Symbol, Any}
end
Bac() = Bac(Dict{Symbol, Any}())

Base.getproperty(x::Bac, property::Symbol) = getfield(x, :properties)[property]
Base.setproperty!(x::Bac, property::Symbol, value) = getfield(x, :properties)[property] = value
Base.propertynames(x::Bac) = keys(getfield(x, :properties))
