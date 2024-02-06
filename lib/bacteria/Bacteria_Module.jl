module Bacteria_Module
"""
This module includese all the files that are in the bacteria map.
The module can be useful when a file requires all of the bacteria files.
"""

using InvertedIndices

include(string(pwd(), "\\lib\\bacteria\\bacteria_detachment.jl"))
include(string(pwd(), "\\lib\\bacteria\\bacteria_die.jl"))
include(string(pwd(), "\\lib\\bacteria\\bacteria_divide.jl"))
include(string(pwd(), "\\lib\\bacteria\\bacteria_inactivate.jl"))
include(string(pwd(), "\\lib\\bacteria\\granule_density.jl"))
include(string(pwd(), "\\lib\\bacteria\\killBacs.jl"))
include(string(pwd(), "\\lib\\bacteria\\sort_bacteria_into_chunks.jl"))
include(string(pwd(), "\\lib\\bacteria\\update_bacterial_mass.jl"))
include(string(pwd(), "\\lib\\bacteria\\update_bacterial_radius.jl"))

end
