function increase_dT_diffusion!(Time, msg, dx, constants)
    """
    This function increases the dT for diffusion timestepping
    
    Arguments
    Time:               A "General" struct containing all time-related parameters
    msg:                The reason for increasing the dT (diffusion)
    dx:                 The size of a gridcell
    constants:          A "General" struct containing all the simulation constants

    Returns
    Time                A Time struct with increased dT (diffusion)
    """
    Time.dT = minimum([Time.dT / 0.9, Time.maxDT])
    Time.changed_dT = Time.current
    println("$(msg), \n \tthus dT increased to $(Time.dT)")
    println("Neumann value of stability: $(maximum(constants.diffusion_rates * Time.dT / (dx ^ 2)))") 

    return Time
end
