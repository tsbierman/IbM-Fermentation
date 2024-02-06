module Lib_Module
using DSP

"""
This module includes most of the files in the lib directory.
This can be useful when a lot of these functions have to be used in the same file
"""

include(string(pwd(), "\\lib\\Struct_Module.jl"))
include(string(pwd(), "\\lib\\calculate_bulk_concentrations.jl"))
include(string(pwd(), "\\lib\\create_chunks.jl"))
include(string(pwd(), "\\lib\\determine_diffusion_region.jl"))
include(string(pwd(), "\\lib\\determine_where_bacteria_in_grid.jl"))
include(string(pwd(), "\\lib\\set_concentrations.jl"))
end
