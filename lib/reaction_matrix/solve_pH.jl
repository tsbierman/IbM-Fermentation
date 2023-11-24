function calculate_spcM!(spcM, Sh, Keq, StV)
    """
    This function calculates the specieMatrix (spcM) based on the compound concentrations,
    the Keq and the proton concentration.
    """

    Denm = (1 .+ Keq[:, 1]) * Sh^3 .+ Keq[:, 2] * Sh^2 .+ Keq[:, 2] .* Keq[:, 3] * Sh .+ Keq[:, 2] .* Keq[:, 3] .* Keq[:, 4] # Common denominator for all equations

    # Calculate concentrations for all the species. Even when not all species are occupied for some compounds,
    # thus some Keq are 0, the general formula will still work. Unnecessary parts will cancel outflow
    spcM[:,1] = (Keq[:,1]) .* StV .* Sh^3 ./ Denm                     # Gaseous equivalent
    spcM[:,2] = (StV .* Sh^3) ./ Denm                                 # Fully protonated
    spcM[:,3] = (StV .* Sh^2 .* Keq[:,2]) ./ Denm                     # First deprotonation
    spcM[:,4] = (StV .* Sh .* Keq[:,2] .* Keq[:,3]) ./ Denm           # Second deprotonation
    spcM[:,5] = (StV .* Keq[:,2] .* Keq[:,3] .* Keq[:,4]) ./ Denm     # Third deprotonation

    return spcM
end


function solve_pH(Sh_ini, StV, Keq, chrM, calculate_pH, Tol)
    """
    Solve the pH and speciation per grid cell using a modified Newton-Raphson algorithm

    Sh_ini: initial guess for the H+ concentration
    StV: concentration vector
    Keq: equilibrium constants per subcompound of the StV vector
    chrM: charge per subcompound of the StV vector
    calculate_pH: (bool) if false, then assume pH is a set value
    Tol: Newton-Raphson tolerance

    Returns:
    spcM: species matrix for the respective pH
    Sh: proton concentration
    """

    if !calculate_pH
        # Assume steady pH value, only calculate speciation
        spcM = zeros(size(chrM))
        Sh = Sh_ini
        spcM = calculate_spcM!(spcM, Sh, Keq, StV) # Fill spcM

    else
        # Use Newton-Raphson method
        Sh = Sh_ini # Initial guess
        ipH = 1 # counter of convergences
        maxIter = 20
        err = 1 # initial error
        F = 1 # Initial error of charge balance

        # Initialisation of matrix of species
        spcM = zeros(size(chrM))
        dspcM = zeros(size(chrM))

        while (abs(err) > Tol) && (abs(F) > Tol) && ipH <= maxIter
            spcM = calculate_spcM!(spcM, Sh, Keq, StV) # Fill spcM based on current Sh

            # Evaluate charge balance for current Sh, F(Sh)
            F = Sh + sum(spcM .* chrM) # Charge balance

            # Calculation of all derivated functions
            Denm = (1 .+ Keq[:, 1]) * Sh^3 .+ Keq[:, 2] * Sh^2 .+ Keq[:, 2] .* Keq[:, 3] * Sh .+ Keq[:, 2] .* Keq[:, 3] .* Keq[:, 4] # Common denominator for all equations
            dDenm = Denm .^2
            aux = 3 * Sh^2 * (Keq[:, 1] .+ 1) .+ 2 * Sh * Keq[:, 2] .+ Keq[:, 2] .* Keq[:, 3]

            dspcM[:,1] = (3 * Sh^2 * Keq[:, 1] .* StV) ./ (Denm) .- ((Keq[:, 1] .* StV * Sh^3) .* aux) ./ (dDenm)
            dspcM[:,2] = (3 * Sh^2 * StV) ./ Denm .- (StV * Sh^3 .* aux) ./ dDenm
            dspcM[:,3] = (2 * Sh * Keq[:, 2] .* StV) ./ Denm .- ((Keq[:, 2] .* StV * Sh^2) .* aux) ./ dDenm
            dspcM[:,4] = (Keq[:, 2] .* Keq[:, 3] .* StV) ./ Denm .- ((Keq[:, 2] .* Keq[:, 3] .* StV * Sh) .* aux) ./ dDenm
            dspcM[:,5] = -(Keq[:, 2] .* Keq[:, 3] .* Keq[:, 4] .* StV .* aux) ./ dDenm

            # Evalutaion of the charge balance for the current Sh value, dF(Sh)
            dF = 1 + sum(spcM .* chrM)
            # Error
            err = F / dF
            # Newton-Raphson
            Sh = maximum([Sh - err, 1e-10])

            ipH = ipH + 1
        end

        if any(spcM .< 0)
            @warn("DEBUG:actionRequired, debug: negative concentration encountered after pH calculation...")
        end

        if (Sh < 1e14 || Sh > 1)
            @warn("DEBUG:actionRequired, debug: pH found outside of 1-14 range...")
        end

        if ipH >= maxIter
            throw(ErrorException("pH solver did not converge to a solution within given number of iterations"))
        end
    end
    return spcM, Sh
    
end
