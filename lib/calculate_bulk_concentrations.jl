function massbal(bulk_conc, p, t)
    """'
    Differential Equation for the mass balance over the entire reactor
    Will modify the HRT to match the setpoint of NH# if the outflow concentration is larger
    than the setpoint.

    Arguments of function is not like other languages, which have (t,y,parameters...)
    This uses the (u,p,t), with u the variables that needs to change over time and p the parameters
    p thus contains a lot of parameters:
        cumulative_reacted: cumulative reaction rate [mol/L/h] per compound
        reactor_influx: influx concentrations [mol/L] per compound
        variableHRT: Boolean indicating whether HRT is variable
        bulk_setpoint: setpoint for outflux of compound[setpoint_index]
        Dir_k: BitArray of which compounds follow Dirichlet boundary condition
        settings: a struct containing all the reactor settings
        invHRT: 1/hydrolic retention time

    returns: dy --> derivative of bulk concentration
    """
    cumulative_reacted = p[1]
    reactor_influx = p[2]
    variableHRT = p[3]
    bulk_setpoint = p[4]
    setpoint_index = p[5]
    Dir_k = p[6]
    settings = p[7]
    invHRT = p[8]

    dy = zeros(length(bulk_conc))

    if settings.structure_model
        # These cases are quite hardcoded, therefore, they are just copied to
        # obtain something complete. I am not sure whether they make sense.
        if settings.type == "Neut"
            # For index 1, if HRT is changeable and bulk concentration is not the set_point
            if variableHRT == true && (bulk_conc[1] > bulk_setpoint || bulk_conc[1] < bulk_setpoint)
                # Only able to adjust when reaction < 0. If reaction > 0, this would mean that Cin - c < 0
                # invHRT can than not be adjusted to fix change at 0, so dy is updated
                if cumulative_reacted[1] < 0
                    invHRT = -cumulative_reacted[1] / (reactor_influx[1] - bulk_setpoint)
                end
            else
                dy[1] = invHRT * (reactor_influx[1] - bulk_conc[1]) + cumulative_reacted[1]
            end
            # For index 2
            if variableHRT == true && (bulk_conc[2] > bulk_setpoint || bulk_conc[2] < bulk_setpoint)
                if cumulative_reacted[2] < 0
                    invHRT = -cumulative_reacted[2] / (reactor_influx[2] - bulk_setpoint)
                end
            else
                dy[2] = invHRT * (reactor_influx[2] - bulk_conc[2]) + cumulative_reacted[2]
            end
            # For index 3
            if variableHRT == true && (bulk_conc[3] > bulk_setpoint || bulk_conc[3] < bulk_setpoint)
                if cumulative_reacted[3] < 0
                    invHRT = -cumulative_reacted[3] / (reactor_influx[3] - bulk_setpoint)
                end
            else
                dy[3] = invHRT * (reactor_influx[3] - bulk_conc[3]) + cumulative_reacted[3]
            end

            dy[4:end] = invHRT * (reactor_influx[4:end] .- bulk_conc[4:end]) .+ cumulative_reacted[4:end]
            dy[5:end] = 0

        elseif settings.type in ("Comp", "Comm", "Copr")
            # For index 1
            if variableHRT == true && (bulk_conc[1] > bulk_setpoint || bulk_conc[1] < bulk_setpoint)
                if cumulative_reacted[1] < 0
                    invHRT = -cumulative_reacted[1] / (reactor_influx[1] - bulk_setpoint)
                end
            else
                dy[1] = invHRT * (reactor_influx[1] - bulk_conc[1]) + cumulative_reacted[1]
            end

            dy[2:end] = invHRT * (reactor_influx[2:end] .- bulk_conc[2:end]) .+ cumulative_reacted[2:end]
            dy[5:end] = 0

        else
            throw(ErrorException("Type <$(settings.type)> is not a registered set of simulations"))
        end

    else
        if variableHRT
            # If bulk_conc < setpoint, decrease HRT to value that would suffice
            if cumulative_reacted[setpoint_index] < 0
                invHRT = -cumulative_reacted[setpoint_index] / (reactor_influx[setpoint_index] - bulk_setpoint)
            end
        end

        # Always calculate the change of non-dirichlet bulk concentrations
        dy[.!Dir_k] = invHRT * (reactor_influx[.!Dir_k] - bulk_conc[.!Dir_k]) .+ cumulative_reacted[.!Dir_k]

    end

    return dy
end


function correct_negative_concentrations(conc)
    """
    Perform a correction to get rid of negative concentrations
    conc: bulk concentrations before correcting
    
    ------------------------------------------------------------
            CONTAINS HARDCODED LINES FOR NITROSPIRA 
    ------------------------------------------------------------

    Returns
    conc: bulk concentrations after correction
    """

    negative_indices = findall(conc .< 0)

    if length(negative_indices) > 0
        warn("DEBUG: noActionRequired, debug: negative concentration encountered and corrected")
        for neg_index in negative_indices
            # if [NH3] < 0, then remove excess consumption from [NO2]
            if neg_index == 1
                conc[2] = conc[2] + conc[1]
                conc[1] = 0

            # if [NO2] < 0, then remove excess consumption from [NO3]
            elseif neg_index == 2
                conc[3] = conc[3] + conc[2]
                conc[2] = 0
            end
        end
    end
end


function controlpH(Keq, chrM, compoundNames, pH, conc)
    """
    Add NaHCO3 to control the pH at a certain point

    Keq: Matrix with Equilibrium constants for all compounds and deprotonations
    chrM: Matrix with charges for all compounds and deprotonations
    compoundNames: Compound Names (without H2O and H)
    pH: setpoint of the pH
    conc: bulk concentration per compound before NaHCO3 addition
    
    Returns:
    Conc: bulk concentration per compound after pH control
    """
    
    Tol = 1e-14         # Very close to 0, alligns with pH 14
    Tp = 1              # Initialisation
    u = [conc; 1; 0]    # Add numbers for H2O and H concentrations (even though H is always empty)
    # u contains the total concentration of the compunds (summation of all species)
    NaHCO3 = conc[findall(compoundNames .== "Na")][1] # Assume all Na comes from NaHCO3
    w = 1                                          # Correction 

    spcM = zeros(size(chrM)) # Store the calculated species
    Sh = 10 ^(-pH)           # Concentration of protons

    while abs(Tp) > Tol
        u[findall(compoundNames .== "Na")] .= NaHCO3         # "Add" NaHCO3 from previous iteration to the system
        u[findall(compoundNames .== "CO2")] .= NaHCO3        # "Add" NaHCO3 from previous iteration to the system

        Denm = (1 .+ Keq[:, 1] / w) * Sh^3 .+ Keq[:, 2] * Sh^2 .+ Keq[:, 2] .* Keq[:, 3] * Sh .+ Keq[:, 2] .* Keq[:, 3] .* Keq[:, 4] # Common denominator for all equations

        # Calculate concentrations for all the species. Even when not all species are occupied for some compounds,
        # thus some Keq are 0, the general formula will still work. Unnecessary parts will cancel outflow
        spcM[:,1] = (Keq[:,1] / w) .* u .* Sh^3 ./ Denm                 # Gaseous equivalent
        spcM[:,2] = (u .* Sh^3) ./ Denm                                 # Fully protonated
        spcM[:,3] = (u .* Sh^2 .* Keq[:,2]) ./ Denm                     # First deprotonation
        spcM[:,4] = (u .* Sh .* Keq[:,2] .* Keq[:,3]) ./ Denm           # Second deprotonation
        spcM[:,5] = (u .* Keq[:,2] .* Keq[:,3] .* Keq[:,4]) ./ Denm     # Third deprotonation
        Tp = Sh + sum(spcM .* chrM)         # Charge balance (protons + charge for all other compounds)
        # Tp will be negative when we are not at the solution yet, thus increasing NaHCO3, mimicking an addition
        NaHCO3 = NaHCO3 - Tp
    end

    # Set Na and CO2 concentrations again for bulk
    conc[findall(compoundNames .== "Na")] .= NaHCO3
    conc[findall(compoundNames .== "CO2")] .= NaHCO3

    return conc
end


function calculate_bulk_concentrations(bac, constants, prev_conc, invHRT, reactionMatrix, dT, settings)
    """
    Function to calculate the bulk layer concentrations. Assumes that the simulated
    bio-aggregate is representative of the entire reactor.

    Includes pH correction using NaHCO3

    Has a correction for when accumulated reation in the aggregate exceeds the available
    supply/concentration of compounds.

    constants is a struct containing all constants
    prev_conc: vector with the previous bulk concentration per compound
    invHRT: 1/ (hydrolic retention time)
    reactionMatrix: Matrix with per grid cell and per compound the change [h-1] due to
    bacterial activity
    dT time of integration with diffusion

    Returns:
    bulk_concentrations: vector with the bulk concentration per compound
    invHRT: new 1 / (hydrolic retention time) [h-1]
    """

    # For easy use: unpack constants
    Keq = constants.Keq                                 # Matrix with Equilibrium constants for all compounds and deprotonations
    chrM = constants.chrM                               # Matrix with charges for all compounds and deprotonations
    compoundNames = constants.compoundNames            # Compound Names (without H2O and H)
    pH = constants.pHsetpoint                           # pH
    Vr = constants.Vr                                   # [L] Representative volume of reactor that is modelled
    Vg = constants.Vg                                   # [L] Volume of a grid cell
    Dir_k = constants.Dir_k                             # Booleans whether compounds follow Dirichlet boundary condition
    influent = constants.influent_concentrations        # [mol/L] influent concentrations
    variableHRT = settings.variableHRT                  # Boolean

    if variableHRT
        bulk_setpoint = constants.bulk_setpoint
        setpoint_index = constants.setpoint_index
    else
        bulk_setpoint = 0
        setpoint_index = 0
    end

    if isa(reactionMatrix, Number) # If no reactionMatrix is formed (initiation), take previous concentrations
        bulk_concentrations = prev_conc
    else

        # Calculate correction factor between volume of sphere and volume of slice (simulation domain)
        if settings.model_type in ("granule", "mature granule")
            x = bac.x[bac.active]
            y = bac.y[bac.active]
            centre_x = mean(x)
            centre_y = mean(y)
            radius_granule = maximum(sqrt.((x .- centre_x) .^2 + (y .- centre_y) .^2))
            f = 4 * radius_granule / (3 * constants.bac_max_radius * 2)
        else
            f = 1
        end

        # Combination of dropdims and sum with those dimensions results in a vector that contains total change
        # over the whole matrix per compound. This is then adjusted to a single grid cell
        cumulative_reacted = dropdims(sum(reactionMatrix, dims=(1,2)), dims=(1,2)) * Vg * f / Vr # [mol/L /h]

        try
            parameters = [cumulative_reacted, reactor_influx, variableHRT, bulk_setpoint, setpoint_index, Dir_k, settings, invHRT]
            prob = ODEProblem(massbal, prev_conc, (0.0, dT), parameters)
            sol = solve(prob, Tsit5(), isoutofdomain=(y,p,t)->any(x->x.<0,y), reltol=1e-8, abstol=1e-20)
            bulk_conc_temp = sol.u[end]
            bulk_concentrations = correct_negative_concentrations(bulk_conc_temp) #<E: Negative concentration from mass balance of reactor. />

        catch e # it should just work
            println("Something bad has happened \n $(typeof(e))")
            rethrow(e)
        end

        # Set Dirichlet boundary condition
        bulk_concentrations[Dir_k] = prev_conc[Dir_k]
    end

    # Apply pH correction to bulk_concentrations
    if settings.pHbulkCorrection
        # Hardcoded stuff
        bulk_concentrations[findall(compoundNames .== "SO4")] = bulk_concentrations[findall(compoundNames .== "A")] ./ 2 # A was originaly NH3, but not present in testfile
        bulk_concentrations = controlpH(Keq, chrM, compoundNames, pH, bulk_concentrations)

        if any(bulk_concentrations .< 0)
            warn("DEBUG:actionRequired, debug: negative bulk concentration encountered after pH control... correction required?")
            bulk_concentrations = bulk_concentrations .* (bulk_concentrations .> 0)
        end
    end

    return bulk_concentrations, invHRT
end
