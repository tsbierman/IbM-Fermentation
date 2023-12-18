function slow_convergence(iRES, RESvalues, constants)
    """
    This function detects whether convergence of the diffusion is slow after a certain
    number of diffusion iterations in the steady-state cycle.

    Arguments
    iRES:               The number of steady state checks since last bacterial advancement (dt_bac change)
    RESvalues:          The RESisdual values of this steady state
    constants:          A "General" struct containing all the simulation constants

    Returns
    slow:               A Boolean indicating whether the convergence is low
    """

    include(string(pwd(), "\\lib\\dynamic_dT\\non_convergent.jl"))

    at_cycle_time = mod(iRES, ceil(constants.dynamicDT.nItersCycle/constants.nDiffusion_per_SScheck)) == 0                                                  # only at nCycle iterations
    if at_cycle_time

        no_conv = non_convergent(iRES, RESvalues, constants.dynamicDT.tolerance_no_convergence)                                                             # Better convergence than no-convergence
        direction_of_conv = maximum(RESvalues[:,iRES-1]) - maximum(RESvalues[:,iRES]) > 0 && maximum(RESvalues[:,iRES-2]) - maximum(RESvalues[:,iRES]) > 0  # There has to be at least some convergence
        little_conv = maximum(RESvalues[:,iRES-1]) - maximum(RESvalues[:,iRES]) < 1e-3 && maximum(RESvalues[:,iRES-2]) - maximum(RESvalues[:,iRES]) < 2e-3  # There has to be some convergence (not too much)
        slow = !no_conv && direction_of_conv && little_conv

    else
        slow = false

    end
    return slow
end
