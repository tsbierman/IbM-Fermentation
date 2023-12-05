function increase_dT_bac!(Time, msg)
    """
    Utility function for increasing the dT for bacterial stepping
    """

    Time.dT_bac = minimum([Time.dT_bac / 0.8, Time.maxDT_bac])
    Time.changed_dT_bac = Time.current
    println("$(msg), \n \tthus dT_bac increased to $(Time.dT_bac)")

    return Time
end
