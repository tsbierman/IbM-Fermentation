#Initialise structs that will have to be used later
# Structs need to be declared at top level, so as this
struct General
    properties::Dict{Symbol, Any}
end
General() = General(Dict{Symbol, Any}())

Base.getproperty(x::General, property::Symbol) = getfield(x, :properties)[property]
Base.setproperty!(x::General, property::Symbol, value) = getfield(x, :properties)[property] = value
Base.propertynames(x::General) = keys(getfield(x, :properties))
