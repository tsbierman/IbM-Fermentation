function calculate_spcM!(spcM, Sh, Keq, StV)
    """
    This function calculates the specieMatrix (spcM) based on the compound concentrations,
    the Keq and the proton concentration.

    Arguments
    spcM:           A (ncompounds, 5) matrix containing the current concentrations per specie
    Sh:             The proton concentration
    Keq:            A (ncompounds,) vector containing all equilibrium constants (Kd and Ka)
    StV:            A (ncompounds,) vector containing the total concentration per compound (so sum(species))

    Returns
    spcM:           A (ncompounds, 5) matrix with updated concentrations per specie
    """

    Denm = (1 .+ Keq[:, 1]) .* Sh^3 .+ Keq[:, 2] .* Sh^2 .+ Keq[:, 2] .* Keq[:, 3] .* Sh .+ Keq[:, 2] .* Keq[:, 3] .* Keq[:, 4] # Common denominator for all equations

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
    This function solves the pH and speciation per gridcell using a modified Newton-Raphson algorithm

    Arguments
    Sh_ini:         The initial guess for the proton concentration
    StV:            A (ncompounds,) vector containing the total concentration per compound (so sum(species))
    Keq:            A (ncompounds,) vector containing all equilibrium constants (Kd and Ka)
    chrM:           A (ncompounds, 5) matrix indicating the charge per specie
    calculate_pH:   A Boolean(bool) indicating whether pH needs to be calculated, if false, pH is a set value
    Tol:            The Newton-Raphson tolerance

    Returns
    spcM:           A (ncompounds, 5) matrix containing the current concentrations per specie for the current pH
    Sh:             The updated proton concentration
    """

    if !Bool(calculate_pH)
        # Assume steady pH value, only calculate speciation
        spcM = zeros(size(chrM))
        Sh = Sh_ini
        spcM = calculate_spcM!(spcM, Sh, Keq, StV) # Fill spcM

    else
        # Use Newton-Raphson method
        Sh = Sh_ini                     # Initial guess
        ipH = 1                         # counter of convergences
        maxIter = 20                    # Max iterations
        err = 1                         # initial error
        F = 1                           # Initial error of charge balance

        # Initialisation of matrix of species
        spcM = zeros(size(chrM))
        dspcM = zeros(size(chrM))

        while (abs(err) > Tol) && (abs(F) > Tol) && ipH <= maxIter
            spcM = calculate_spcM!(spcM, Sh, Keq, StV)                      # Fill spcM based on current Sh

            # Evaluate charge balance for current Sh, so evaluate F(Sh)
            F = Sh + sum(spcM .* chrM)                                      # Charge balance

            # Calculation of all derivated functions
            Denm = (1 .+ Keq[:, 1]) .* Sh^3 .+ Keq[:, 2] .* Sh^2 .+ Keq[:, 2] .* Keq[:, 3] .* Sh .+ Keq[:, 2] .* Keq[:, 3] .* Keq[:, 4] # Common denominator for all equations
            dDenm = Denm .^2
            aux = 3 .* Sh^2 .* (Keq[:, 1] .+ 1) .+ 2 .* Sh .* Keq[:, 2] .+ Keq[:, 2] .* Keq[:, 3]

            dspcM[:,1] = (3 .* Sh^2 .* Keq[:, 1] .* StV) ./ (Denm) .- ((Keq[:, 1] .* StV .* Sh^3) .* aux) ./ (dDenm)
            dspcM[:,2] = (3 .* Sh^2 .* StV) ./ Denm .- (StV .* Sh^3 .* aux) ./ dDenm
            dspcM[:,3] = (2 .* Sh .* Keq[:, 2] .* StV) ./ Denm .- ((Keq[:, 2] .* StV .* Sh^2) .* aux) ./ dDenm
            dspcM[:,4] = (Keq[:, 2] .* Keq[:, 3] .* StV) ./ Denm .- ((Keq[:, 2] .* Keq[:, 3] .* StV .* Sh) .* aux) ./ dDenm
            dspcM[:,5] = -(Keq[:, 2] .* Keq[:, 3] .* Keq[:, 4] .* StV .* aux) ./ dDenm

            # Evalutaion of the charge balance for the current Sh value, so evaluate dF(Sh)
            dF = 1 + sum(dspcM .* chrM)
            # Error
            err = F / dF
            # Newton-Raphson
            Sh = maximum([Sh - err, 1e-10])

            ipH = ipH + 1
        end

        # Perform some checks
        if any(spcM .< 0)
            @warn("DEBUG:actionRequired, debug: negative concentration encountered after pH calculation...")
        end
        
        if (Sh < 1e-14 || Sh > 1)
            @warn("DEBUG:actionRequired, debug: pH found outside of 1-14 range...")
        end

        if ipH >= maxIter
            throw(ErrorException("pH solver did not converge to a solution within given number of iterations"))
        end
    end
    return spcM, Sh
    
end
