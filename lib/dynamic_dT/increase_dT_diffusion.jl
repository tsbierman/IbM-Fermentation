function increase_dT_diffusion!(Time, msg, dx, constants_vecfloat)
    """
    This function increases the dT for diffusion timestepping
    
    Arguments
    Time:               A "Float" struct containing time-related parameters of type Float64
    msg:                The reason for increasing the dT (diffusion) (String)
    dx:                 The size of a gridcell
    constants_vecfloat: A "VectorFloat" struct containing simulation constants of type Vector{Float64}

    Returns
    Time:               A "Float" struct containing time-related parameters of type Float64 with increased dT (diffusion)
    """
    
    Time.dT = minimum([Time.dT / 0.9, Time.maxDT])
    Time.changed_dT = Time.current

    @printf("%s, \n \tthus dT increased to %.3e\n", msg, Time.dT)
    @printf("Neumann value of stability: %.5e\n", maximum(constants_vecfloat.diffusion_rates .* Time.dT ./ (dx .^ 2)))
    
    return Time
end
