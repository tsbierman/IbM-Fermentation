function isReached_compound(RES, method, steadystate_tolerance)
    """
    Determine whether steadystate is reached for one compound

    Arguments
    RES:                    A vector that contains the residual of the diffusion-reaction equation for each gridcell
                            in the diffusion region for 1 compound
    method:                 the method which is used to determine steady state
    steadystate_tolerance:  The tolerance when it comes to steadystate

    Returns
    SSreached:              A Boolean whether steady state is reached for this compound
    """

    if method == "mean"
        # Paraphrased from Chiel: When using suspension mode (little granules), this method is not
        # meaningful anymore as the bulk layer has by definition a RES of 0. As the domain in this mode
        # is not shrunk down, a lot of bulk liquid points will be taken, thus lowering the mean. Now
        # instead only the mean of non-zero entries is taken.
        SSdif = mean(abs.(RES[find(RES)])) # Will only return NaN or missing if RES contains NaN or missing

    elseif method == "max"
        SSdif = maximum(abs.(RES))

    elseif method == "norm"
        SSdif = sqrt(sum(RES .^2)) # Calculates as Root-Sum-Squares (RSS)
        # C:in testing this was too harsh of a constraint, because all
        # little errors accumulate and make it so that it takes
        # (figuratively) a million iterations before steady state...

    else
        @error("RES method <$(method)> is not a valid method")

    end
    SSreached = SSdif <= steadystate_tolerance
    return SSreached
end


function steadystate_is_reached(conc, reaction_matrix, dx, bulk_concentrations, diffRegion, constants)
    """
    Checks whether a steady state is reached during the diffusion

    Arguments
    conc:                   A (ny, nx, ncompounds) matrix that includes all concentrations for these gridcells
    reaction_matrix:        A (ny, nx, ncompounds) matrix that contains how much each concentrations changes due to activity
    dx:                     The size of a grid cell
    bulk_concentrations:    A (ncompounds,) vector containing the bulk concentrations
    diffRegion:             A (ny, nx) BitMatrix indicating whether each gridcell is in the diffusion region.
    constants:              A "General" struct containing all the simulation constants

    Returns
    isReached:              A Boolean indicating whether steadystate is reached
    max_RES_value:          A (ncompounds) vector that contains the maximum residual value per compound
    """

    include(string(pwd(), "\\lib\\diffusion\\create_dirichlet_boundary.jl"))

    steadystate_tolerance = constants.steadystate_tolerance     # [mol/L/h]
    method = constants.RESmethod                                # mean, max or norm

    L = [0 1 0; 1 -4 1; 0 1 0]                                  # 2D Laplacian stencil base
    nCompounds = length(constants.compoundNames)
    characteristic_time = dx^2 ./ constants.diffusion_rates     # [h]
    compound_steadystate = BitArray(zeros(nCompounds))          # Binary storage system

    max_RES_value = zeros(nCompounds)

    # per compound, calculate change in concentration
    for index in 1:nCompounds

        # Get extra layer around matrix for correct diffusion at boundary points
        padded_conc = create_dirichlet_boundary(conc[:,:, index], bulk_concentrations[index])   # [mol/L]

        # Use convolution to calculate the diffusion, but only keep the original matrix shape
        RES = conv(padded_conc, L)[3:end-2,3:end-2]                                             # [mol/L]
        RES = RES ./ characteristic_time[index] + reaction_matrix[:,:,index]                    # [mol/L/h]

        # Only look at RES value in the diffusion region, RES of bulk layer has by definition 0 RES.
        compound_steadystate[index] = isReached_compound(RES[diffRegion], method, steadystate_tolerance)
        max_RES_value[index] = maximum(abs.(RES[diffRegion]))                                   # If method = max, same calculation
    end

    isReached = all(compound_steadystate)

    return isReached, max_RES_value
end
