function decrease_dT_bac!(Time, msg)
    """
    Utility function for decreasing the dT for bacterial stepping
    """

    Time.dT_bac = maximum([Time.dT_bac * 0.8, Time.minDT_bac])
    Time.changed_dT_bac = Time.current
    println("$(msg), \n \tthus dT_bac decreased to $(Time.dT_bac)")

    return Time
end
