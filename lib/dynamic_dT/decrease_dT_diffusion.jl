function decrease_dT_diffusion!(Time, msg, dx, constants_vecfloat)
    """
    This function decreases the dT for diffusion timestepping
    
    Arguments
    Time:               A "General" struct containing all time-related parameters
    msg:                The reason for decreasing the dT (diffusion)
    dx:                 The size of a gridcell
    constants:          A "General" struct containing all the simulation constants

    Returns
    Time                A Time struct with decreased dT (diffusion)
    """
    
    prev = Time.dT

    Time.dT = maximum([Time.dT * 0.9, Time.minDT])
    if Time.dT != prev # only set changed_dT if time actually changed_dT
        Time.changed_dT = Time.current
    end
    @printf("%s, \n \tthus dT decreased to %.3e\n", msg, Time.dT)
    @printf("Neumann value of stability: %.5e\n", maximum(constants_vecfloat.diffusion_rates .* Time.dT ./ (dx .^ 2)))

    return Time
end
