function increase_dT_bac!(Time, msg)
    """
    This function increases the dT for bacterial timestepping

    Arguments
    Time:               A "Float" struct containing time-related parameters of type Float64
    msg:                The reason for increasing the dT_bac (String)

    Returns
    Time:               A "Float" struct containing time-related parameters of type Float64 with increased dT_bac
    """

    Time.dT_bac = minimum([Time.dT_bac / 0.8, Time.maxDT_bac])
    Time.changed_dT_bac = Time.current
    @printf("%s, \n \tthus dT_bac increased to %.3e\n", msg, Time.dT_bac)

    return Time
end
