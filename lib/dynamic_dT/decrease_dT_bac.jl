function decrease_dT_bac!(Time, msg)
    """
    This function decreases the dT for bacterial timestepping

    Arguments
    Time:           A "General" struct containing all time-related parameters
    msg:            The reason for decreasing the dT_bac

    Returns
    Time            A Time struct with decreased dT_bac
    """

    Time.dT_bac = maximum([Time.dT_bac * 0.8, Time.minDT_bac])
    Time.changed_dT_bac = Time.current
    println("$(msg), \n \tthus dT_bac decreased to $(Time.dT_bac)")

    return Time
end
