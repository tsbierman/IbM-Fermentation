function upward_trend(iRES, RESvalues)
    """
    This function detects whether there is an upward trend in the RES values
    
    Arguments
    iRES:               The number of steady state checks since last bacterial advancement (dt_bac change)
    RESvalues:          The RESisdual values of this steady state

    Returns
    upward:             A Boolean indicating whether there is an upward trend in thre RES values
    """

    upward = iRES > 2 && any(RESvalues[:, iRES] .- RESvalues[:,iRES-1] .> 0.1) # Checks whether the RES massively diverges
    return upward
end
