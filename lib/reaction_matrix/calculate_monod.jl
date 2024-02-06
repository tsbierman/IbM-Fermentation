function calculate_monod(Ks, Ki, conc)
    """
    This function calculates the Monod coefficient given the Ks's, Ki's and the corresponding concentrations
    
    Arguments
    Ks:         A (ncompounds,) vector with all Ks values
    Ki:         A (ncompounds,) vector with all Ki values
    conc:       A (ncompounds,) vector with all concentrations of the reactive species

    Returns
    M:          A value, which is the Monod part for adjusting the maximal growth
    """

    # Apply Ks
    M = prod( (conc .+ 1e-25) ./ (conc .+ Ks .+ 1e-25))     # + 1e-25 to prevent NaN when conc == 0 and Ks == 0

    # Apply Ki
    for index in eachindex(Ki)
        if Ki[index] != 0
            M = M * Ki[index] / (Ki[index] + conc[index])
        end
    end

    return M
end
