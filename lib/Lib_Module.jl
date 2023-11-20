module Lib_Module
using DSP

include(string(pwd(), "\\lib\\Struct_Module.jl"))

calc_bulk_file = string(pwd(), "\\lib\\calculate_bulk_concentrations.jl")
include(calc_bulk_file)

create_chunks_file = string(pwd(), "\\lib\\create_chunks.jl")
include(create_chunks_file)

deter_diff_region_file = string(pwd(), "\\lib\\determine_diffusion_region.jl")
include(deter_diff_region_file)

deter_bac_grid_file = string(pwd(), "\\lib\\determine_where_bacteria_in_grid.jl")
include(deter_bac_grid_file)

set_conc_file = string(pwd(), "\\lib\\set_concentrations.jl")
include(set_conc_file)

end
