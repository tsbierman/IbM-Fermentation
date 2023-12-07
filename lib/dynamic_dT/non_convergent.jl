function non_convergent(iRES, RESvalues, tol)
    """
    Detect whether the diffusion is no longer converging after a concentration
    number of diffusion iteration in the steady-state cycle

    iRES:           number steady state checks since last bacterial advancement
    RESvalues:      The RESisdual values of this steady state
    tol:            Difference between consecutive maximum RES values below which the system is considered non-convergent
    returns         no_convergence, a boolean
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
