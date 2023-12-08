function upward_trend(iRES, RESvalues)
    """
    Detect whether there is an upward trend in the RES values
    """
    upward = iRES > 2 && any(RESvalues[:, iRES] .- RESvalues[:,iRES-1] .> 0.1) # Checks whether the RES massively diverges
    return upward
end
