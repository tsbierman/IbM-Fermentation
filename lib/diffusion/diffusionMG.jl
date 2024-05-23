function calculate_rhs_dirichlet(phi, L_rhs, value, diffRegion)
    """
    Calculate part of the rhs of the diffusion-reaction equation using convolution on 1 compound only.
    and only on the diffusion region. It assumes that there is at least one layer of bulk
    liquid around the diffusion region.

    Arguments
    phi:            A (ny, nx) matrix with starting values for 1 compound
    L_rhs:          A kernel used for the rhs convolution
    value:          The boundary value for this compound
    diffRegion:     A BitMatrix with per gridcell whether that gridcell is in the diffusion region

    Only required if else the concentration would be 0 (outside of domain)

    Returns
    rhs:            The right-hand-side of the diffusion equation due to diffusion, valid only in the diffusion region.
                    Outside the region, an artifical vlaue of the bulk_concentration is set.
    """

    adjusted_phi = diffRegion .* phi .+ .!diffRegion .* value                            # Fix outside diffusion region
    # rhs_diffRegion = conv(adjusted_phi, L_rhs)[2:end-1,2:end-1]                         # Convolution
    rhs_diffRegion = imfilter(adjusted_phi, reflect(centered(L_rhs)), Fill(0))                         # Convolution
    rhs = diffRegion .* rhs_diffRegion .+ .!diffRegion .* ones(size(phi)) .* value       # Fix outside diffusion region

    return rhs
end


function diffusionMG!(conc, reaction_matrix, bulk_concentrations, diffRegion, grid_float, constants_float, constants_vecfloat, constants_vecstring, Time)
    """
    Solve diffusion for all molecules in the liquid phase using the multigrid method.
    IMPORTANT: only runs for all dirichlet conditions as of now. Future versions should
    include variable conditions per boundary.

    Arguments
    conc:                   A (ny, nx, ncompounds) matrix containing the concentration of each molecule in each gridcell [mol/L] 
    reaction_matrix:        A (ny, nx, ncompounds) matrix with per gridcell and per compound the change [mol/L] due to bacterial activity
    bulk_concentrations:    A (ncompounds,) vector with the bulk concentrations per compound
    diffRegion:             A BitMatrix with per gridcell whether that gridcell is in the diffusion region
    grid:                   A "General" struct containing all parameters related to the grid
    constants:              A "General" struct containing all the simulation constants
    Time:                   A "General" struct containing all time-related variables

    Returns
    concs:                  A (ny, nx, ncompounds) matrix containing the concentrations after solving the diffusion equations [mol/L]
    """

    # variable declarations/unpacking
    diffusion_coef = constants_vecfloat.diffusion_rates # [m2/h]
    accuracy = constants_float.diffusion_accuracy
    nCompounds = length(diffusion_coef)                 # Only liquid compounds
    dx = grid_float.dx
    dT = Time.dT

    # Set parameters for V_cycle
    iter_pre = 6
    iter_post = 7
    iter_final = 5

    #convert concentration to mol/m3
    conc = conc .* 1000

    # Stencil Initialisation
    L = [0 1 0; 1 -4 1; 0 1 0]                  # Laplacian stencil
    I = zeros(3,3)
    I[2,2] = 1

    base = [1 2 1; 2 4 2; 1 2 1]
    L_restriction = base ./ 16
    L_prolongation = base ./ 4

    for icompound in 1:nCompounds
        # stencil updates/declarations
        alpha = dT * diffusion_coef[icompound] / (2 * dx^2)
        L_lhs = I .- alpha .* L                  # Lefthand-side stencil
        L_0 = alpha .* L                         # Basis Laplacian stencil
        L_rhs = I .+ alpha .* L                  # Righthand-side stencil

        # Create right hand side
        # - boundary conditions
        rhs_bc = calculate_rhs_dirichlet(conc[:,:,icompound], L_rhs, bulk_concentrations[icompound] .* 1000, diffRegion)     # [ny,nx,ncompounds]

        # - Reaction matrix
        rhs_react = dT .* 1000 .* reaction_matrix[:,:,icompound]

        # - Combine
        rhs = rhs_bc .+ rhs_react

        # solve using multigrid
        isSolution = false      # without running, not isSolution
        while !isSolution
            conc[:,:,icompound] = V_Cycle!(conc[:,:,icompound], diffRegion, bulk_concentrations[icompound] .* 1000, rhs, L_0, L_restriction, L_prolongation, 9, 0, iter_pre, iter_post, iter_final)
            residual_diffRegion = residual(conc[:,:,icompound], rhs, L_lhs)     # Check residuals
            residual_diffRegion = diffRegion .* residual_diffRegion             # Only for the diffusion region
            isSolution = sum(residual_diffRegion .^2) < accuracy^2              # Check whether condition is met, not the case? --> Another V-Cycle
        end

        # CHeck for negative values, raise error
        negative_concentration = conc[:,:,icompound] .< 0
        if any(negative_concentration)
            if Time.dT != Time.minDT && abs(minimum(conc[negative_concentration, icompound])) > accuracy^2 / 100   # If dT can be reduced and is a significant negative value
                throw(ErrorException("Diffusion:negative_concentration, Negative concentration encountered in diffusion solution of compound $(constants_vecstring.compoundNames[icompound])"))
            else                                                                                        # dT cannot be reduced, thus return corrected concentration or insignificant negative value
                temp = conc[:,:,icompound]
                @warn("Diffusion:negative_concentration, Negative concentration encountered in diffusion solution of compound $(constants_vecstring.compoundNames[icompound]), but cannot correct dT value thus corrected $(sum(negative_concentration)) value(s) (smallest number $(minimum(temp[negative_concentration]))) to 0)")
                conc[:,:,icompound] = (conc[:,:,icompound] .> 0) .* conc[:,:,icompound]
            end
        end
    end

    # convert concentration back to mol/L
    conc = conc ./ 1000
   
    return conc
end
