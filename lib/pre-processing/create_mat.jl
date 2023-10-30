#Initialise structs that will have to be used later
# Structs need to be declared at top level

loading_file = string(Base.source_dir(), "\\","loadPresetFile.jl")
include(loading_file)


struct General
    properties::Dict{Symbol, Any}
end
General() = General(Dict{Symbol, Any}())

Base.getproperty(x::General, property::Symbol) = getfield(x, :properties)[property]
Base.setproperty!(x::General, property::Symbol, value) = getfield(x, :properties)[property] = value
Base.propertynames(x::General) = keys(getfield(x, :properties))

import XLSX

# This needs to be changed if the test_file is moved to another directory.
code_folder = dirname(dirname(Base.source_dir()))
filename = string(code_folder, "\\planning\\test_file.xlsx")

grid, bac_init, constants, settings, init_params = loadPresetFile(filename)
