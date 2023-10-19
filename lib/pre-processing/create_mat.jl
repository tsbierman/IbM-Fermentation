#Initialise structs that will have to be used later
# Structs need to be declared at top level, so as this
struct General
    properties::Dict{Symbol, Any}
end
General() = General(Dict{Symbol, Any}())

Base.getproperty(x::General, property::Symbol) = getfield(x, :properties)[property]
Base.setproperty!(x::General, property::Symbol, value) = getfield(x, :properties)[property] = value
Base.propertynames(x::General) = keys(getfield(x, :properties))

import XLSX
filename = string(Base.source_dir(), "\\","test_file.xlsx")
# Should later look something like this:
# filename = ".\\Documents\\Universiteit\\Master\\Jaar 2\\MEP\\Code\\IbM-Fermentation\\IbM-Fermentation\\lib\\pre-processing\\test_file.xlsx"
