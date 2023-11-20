module Bacteria_Module

detachment_file = string(pwd(), "\\lib\\bacteria\\bacteria_detachment.jl")
include(detachment_file)

die_file = string(pwd(), "\\lib\\bacteria\\bacteria_die.jl")
include(die_file)

divide_file = string(pwd(), "\\lib\\bacteria\\bacteria_divide.jl")
include(divide_file)

inactivate_file = string(pwd(), "\\lib\\bacteria\\bacteria_inactivate.jl")
include(inactivate_file)

granule_density_file = string(pwd(), "\\lib\\bacteria\\granule_density.jl")
include(granule_density_file)

using InvertedIndices
killBacs_file = string(pwd(), "\\lib\\bacteria\\killBacs.jl")
include(killBacs_file)

sort_bacteria_file = string(pwd(), "\\lib\\bacteria\\sort_bacteria_into_chunks.jl")
include(sort_bacteria_file)

update_mass_file = string(pwd(), "\\lib\\bacteria\\update_bacterial_mass.jl")
include(update_mass_file)

update_radius_file = string(pwd(), "\\lib\\bacteria\\update_bacterial_radius.jl")
include(update_radius_file)

end
