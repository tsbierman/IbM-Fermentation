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

struct VectorFloat_struct
    properties::Dict{Symbol, Vector{Float64}}
end
VectorFloat_struct() = VectorFloat_struct(Dict{Symbol, Vector{Float64}}())

Base.getproperty(x::VectorFloat_struct, property::Symbol) = getfield(x, :properties)[property]
Base.setproperty!(x::VectorFloat_struct, property::Symbol, value) = getfield(x, :properties)[property] = value
Base.propertynames(x::VectorFloat_struct) = keys(getfield(x, :properties))

struct MatrixFloat_struct
    properties::Dict{Symbol, Matrix{Float64}}
end
MatrixFloat_struct() = MatrixFloat_struct(Dict{Symbol, Matrix{Float64}}())

Base.getproperty(x::MatrixFloat_struct, property::Symbol) = getfield(x, :properties)[property]
Base.setproperty!(x::MatrixFloat_struct, property::Symbol, value) = getfield(x, :properties)[property] = value
Base.propertynames(x::MatrixFloat_struct) = keys(getfield(x, :properties))

struct Int_struct
    properties::Dict{Symbol, Int32}
end
Int_struct() = Int_struct(Dict{Symbol, Int32}())

Base.getproperty(x::Int_struct, property::Symbol) = getfield(x, :properties)[property]
Base.setproperty!(x::Int_struct, property::Symbol, value) = getfield(x, :properties)[property] = value
Base.propertynames(x::Int_struct) = keys(getfield(x, :properties))

struct VectorInt_struct
    properties::Dict{Symbol, Vector{Int32}}
end
VectorInt_struct() = VectorInt_struct(Dict{Symbol, Vector{Int32}}())

Base.getproperty(x::VectorInt_struct, property::Symbol) = getfield(x, :properties)[property]
Base.setproperty!(x::VectorInt_struct, property::Symbol, value) = getfield(x, :properties)[property] = value
Base.propertynames(x::VectorInt_struct) = keys(getfield(x, :properties))

struct MatrixInt_struct
    properties::Dict{Symbol, Matrix{Int32}}
end
MatrixInt_struct() = MatrixInt_struct(Dict{Symbol, Matrix{Int32}}())

Base.getproperty(x::MatrixInt_struct, property::Symbol) = getfield(x, :properties)[property]
Base.setproperty!(x::MatrixInt_struct, property::Symbol, value) = getfield(x, :properties)[property] = value
Base.propertynames(x::MatrixInt_struct) = keys(getfield(x, :properties))

struct Bool_struct
    properties::Dict{Symbol, Bool}
end
Bool_struct() = Bool_struct(Dict{Symbol, Bool}())

Base.getproperty(x::Bool_struct, property::Symbol) = getfield(x, :properties)[property]
Base.setproperty!(x::Bool_struct, property::Symbol, value) = getfield(x, :properties)[property] = value
Base.propertynames(x::Bool_struct) = keys(getfield(x, :properties))

struct VectorBool_struct
    properties::Dict{Symbol, Vector{Bool}}
end
VectorBool_struct() = VectorBool_struct(Dict{Symbol, Vector{Bool}}())

Base.getproperty(x::VectorBool_struct, property::Symbol) = getfield(x, :properties)[property]
Base.setproperty!(x::VectorBool_struct, property::Symbol, value) = getfield(x, :properties)[property] = value
Base.propertynames(x::VectorBool_struct) = keys(getfield(x, :properties))

struct MatrixBool_struct
    properties::Dict{Symbol, Matrix{Bool}}
end
MatrixBool_struct() = MatrixBool_struct(Dict{Symbol, Matrix{Bool}}())

Base.getproperty(x::MatrixBool_struct, property::Symbol) = getfield(x, :properties)[property]
Base.setproperty!(x::MatrixBool_struct, property::Symbol, value) = getfield(x, :properties)[property] = value
Base.propertynames(x::MatrixBool_struct) = keys(getfield(x, :properties))

struct String_struct
    properties::Dict{Symbol, String}
end
String_struct() = String_struct(Dict{Symbol, String}())

Base.getproperty(x::String_struct, property::Symbol) = getfield(x, :properties)[property]
Base.setproperty!(x::String_struct, property::Symbol, value) = getfield(x, :properties)[property] = value
Base.propertynames(x::String_struct) = keys(getfield(x, :properties))

struct VectorString_struct
    properties::Dict{Symbol, Vector{String}}
end
VectorString_struct() = VectorString_struct(Dict{Symbol, Vector{String}}())

Base.getproperty(x::VectorString_struct, property::Symbol) = getfield(x, :properties)[property]
Base.setproperty!(x::VectorString_struct, property::Symbol, value) = getfield(x, :properties)[property] = value
Base.propertynames(x::VectorString_struct) = keys(getfield(x, :properties))
