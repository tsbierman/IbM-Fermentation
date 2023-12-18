function non_convergent(iRES, RESvalues, tol)
    """
    This function detects whether the diffusion is no longer converging after a
    certain number of diffusion iterations in the steady-state cycle

    Arguments
    iRES:               The number of steady state checks since last bacterial advancement (dt_bac change)
    RESvalues:          The RESisdual values of this steady state
    tol:                The difference between consecutive maximum RES values below which the system is considered non-convergent

    Returns
    no_convergence:     A Boolean indicating whether the diffusion is no longer converging
    """

    # Currently, a SS check is done every 2 diffusion iterations
    if iRES > 5 # If lower than 5, not enough steady state checks performed to consider non-convergent
        no_conv_positive = abs(maximum(RESvalues[:, iRES]) - maximum(RESvalues[:, iRES-1])) < tol &&    # Checks difference between current and previous maximum(RES)
        abs(maximum(RESvalues[:, iRES]) - maximum(RESvalues[:, iRES-2])) < 2*tol                        # Checks difference between current and pre-previous maximum(RES), together they check whether convergence goes fast enough

        no_conv_negative = maximum(RESvalues[:, iRES]) - maximum(RESvalues[:, iRES-1]) > 0 &&           # Checks whether current error (iRES) is larger than previous
        abs(maximum(RESvalues[:, iRES]) - maximum(RESvalues[:, iRES-1])) < 1e-5                         # Checks difference between current and previous maximum(RES), together they check whether it diverges

        no_convergence = no_conv_positive || no_conv_negative
        
    else
        no_convergence = false
        
    end
    return no_convergence
end
