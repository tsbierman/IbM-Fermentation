function increase_dT_diffusion!(Time, msg, dx, constants)
    """
    Utility function for increasing the dT for diffusion.
    """

    Time.dT = minimum([Time.dT / 0.9, Time.maxDT])
    Time.changed_dT = Time.current
    println("$(msg), \n \tthus dT increased to $(Time.dT)")
    println("Neumann value of stability: $(maximum(constants.diffusion_rates * Time.dT / (dx ^ 2)))") 

    return Time
end
