# pre_processing files
include("lib\\pre_processing\\create_mat.jl")
include("lib\\pre_processing\\loadPresetFile.jl")
println(">>>>>>>>>>>>>>INCLUDED PRE_PROCESSING FILES")

# General files
include("IbM.jl")
include("extract_data.jl")
println(">>>>>>>>>>>>>>INCLUDED IbM")

# Lib files
include("lib\\calculate_bulk_concentrations.jl")
include("lib\\create_chunks.jl")
include("lib\\determine_diffusion_region.jl")
include("lib\\determine_where_bacteria_in_grid.jl")
include("lib\\initTime.jl")
include("lib\\integTime.jl")
include("lib\\set_concentrations.jl")
include("lib\\Struct_Module.jl")
println(">>>>>>>>>>>>>>INCLUDED LIB FILES")

# Bacteria files
include("lib\\bacteria\\bacteria_detachment.jl")
include("lib\\bacteria\\bacteria_die.jl")
include("lib\\bacteria\\bacteria_divide.jl")
include("lib\\bacteria\\bacteria_inactivate.jl")
include("lib\\bacteria\\bacteria_shove.jl")
include("lib\\bacteria\\granule_density.jl")
include("lib\\bacteria\\killBacs.jl")
include("lib\\bacteria\\sort_bacteria_into_chunks.jl")
include("lib\\bacteria\\update_bacterial_mass.jl")
include("lib\\bacteria\\update_bacterial_radius.jl")
println(">>>>>>>>>>>>>>INCLUDED BACTERIA FILES")

# Detachment files
include("lib\\detachment\\calcTimeOfDetach.jl")
include("lib\\detachment\\calculateLocalDetachmentRate.jl")
include("lib\\detachment\\computeRoot.jl")
include("lib\\detachment\\recalculateT.jl")
println(">>>>>>>>>>>>>>INCLUDED DETACHMENT FILES")

# Diffusion files
include("lib\\diffusion\\create_dirichlet_boundary.jl")
include("lib\\diffusion\\diffusionMG.jl")
include("lib\\diffusion\\prolongation.jl")
include("lib\\diffusion\\residual.jl")
include("lib\\diffusion\\restriction.jl")
include("lib\\diffusion\\smoothing.jl")
include("lib\\diffusion\\steadystate_is_reached.jl")
include("lib\\diffusion\\V_cycle.jl")
println(">>>>>>>>>>>>>>INCLUDED DIFFUSION FILES")

# Dynamic_dT files
include("lib\\dynamic_dT\\bulk_conc_diff_within_limit.jl")
include("lib\\dynamic_dT\\decrease_dT_bac.jl")
include("lib\\dynamic_dT\\decrease_dT_diffusion.jl")
include("lib\\dynamic_dT\\increase_dT_bac.jl")
include("lib\\dynamic_dT\\increase_dT_diffusion.jl")
include("lib\\dynamic_dT\\multiple_high_iters.jl")
include("lib\\dynamic_dT\\multiple_low_initRES.jl")
include("lib\\dynamic_dT\\non_convergent.jl")
include("lib\\dynamic_dT\\slow_convergence.jl")
include("lib\\dynamic_dT\\upward_trend.jl")
println(">>>>>>>>>>>>>>INCLUDED DYNAMIC_DT FILES")

# Post_processing files
include("lib\\post_processing\\save_backup.jl")
include("lib\\post_processing\\save_profile.jl")
include("lib\\post_processing\\save_profiling.jl")
include("lib\\post_processing\\save_slice.jl")
println(">>>>>>>>>>>>>>INCLUDED POST_PROCESSING FILES")

# Reaction_matrix files
include("lib\\reaction_matrix\\calculate_monod.jl")
include("lib\\reaction_matrix\\calculate_reaction_matrix.jl")
include("lib\\reaction_matrix\\determine_max_growth_rate_and_maint.jl")
include("lib\\reaction_matrix\\rMatrix_section.jl")
include("lib\\reaction_matrix\\solve_pH.jl")
println(">>>>>>>>>>>>>>INCLUDED REACTION_MATRIX FILES")

# For testing
include("test\\plotLogicalGrid.jl") # For visual testing

println("DONE INCLUDING")
