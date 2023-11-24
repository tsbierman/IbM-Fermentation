function calculate_monod(Ks, Ki, conc)
    """
    Calculate the Monod coefficient given the Ks's, Ki's and the corresponding concentrations
    
    Ks: vector with all Ks values (1-by-n)
    Ki: vector with all Ki values (1-by-n)
    conc: vector with all concentrations corresponding to the Ks and Ki values (1-by-n)

    M: Monod part for applying towards growth
    """
    # Apply Ks
    M = prod( (conc .+ 1e-25) ./ (conc .+ Ks .+ 1e-25)) # + 1e-25 to prevent NaN when conc == 0 and Ks == 0

    # Apply Ki
    for index in eachindex(Ki)
        if Ki[index] != 0
            M = M * Ki[index] / (Ki[idnex] + conc[index])
        end
    end

    return M
end
