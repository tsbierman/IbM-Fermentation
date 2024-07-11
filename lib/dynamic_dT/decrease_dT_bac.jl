function decrease_dT_bac!(Time, msg)
    """
    This function decreases the dT for bacterial timestepping

    Arguments
    Time:               A "Float" struct containing time-related parameters of type Float64
    msg:                The reason for decreasing the dT_bac (String)

    Returns
    Time:               A "Float" struct containing time-related parameters of type Float64 with decreased dT_bac
    """

    Time.dT_bac = maximum([Time.dT_bac * 0.8, Time.minDT_bac])
    Time.changed_dT_bac = Time.current
    @printf("%s, \n \tthus dT_bac decreased to %.3e\n", msg, Time.dT_bac)

    return Time
end
