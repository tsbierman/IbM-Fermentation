"""
This file contains the definition of the General struct
This struct is used to store variables and be able to 
dynamically add or change its properties 
"""

struct General
    properties::Dict{Symbol, Any}
end
General() = General(Dict{Symbol, Any}())

Base.getproperty(x::General, property::Symbol) = getfield(x, :properties)[property]
Base.setproperty!(x::General, property::Symbol, value) = getfield(x, :properties)[property] = value
Base.propertynames(x::General) = keys(getfield(x, :properties))

struct Float_struct
    properties::Dict{Symbol, Float64}
end
Float_struct() = Float_struct(Dict{Symbol, Float64}())

Base.getproperty(x::Float_struct, property::Symbol) = getfield(x, :properties)[property]
Base.setproperty!(x::Float_struct, property::Symbol, value) = getfield(x, :properties)[property] = value
Base.propertynames(x::Float_struct) = keys(getfield(x, :properties))

struct Int_struct
    properties::Dict{Symbol, Int32}
end
Int_struct() = Int_struct(Dict{Symbol, Int32}())

Base.getproperty(x::Int_struct, property::Symbol) = getfield(x, :properties)[property]
Base.setproperty!(x::Int_struct, property::Symbol, value) = getfield(x, :properties)[property] = value
Base.propertynames(x::Int_struct) = keys(getfield(x, :properties))
