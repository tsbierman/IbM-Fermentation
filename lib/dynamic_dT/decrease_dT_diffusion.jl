function decrease_dT_diffusion!(Time, msg, dx, constants)
    """
    Utility function for decreasing the dT for diffusion
    """
    
    prev = Time.dT

    Time.dT = maximum(Time.dT * 0.9, Time.minDT)
    if Time.dT != prev # only set changed_dT if time actually changed_dT
        Time.changed_dT = Time.current
    end
    println("$(msg), \n \tthus dT decreased to $(Time.dT)")
    println("Neumann value of stability: $(maximum(constants.diffusion_rates * Time.dT / (dx ^ 2)))") 
    return Time
end
