function V_Cycle!(phi, diffRegion, bulk_value, f, L_0, L_restriction, L_prolongation, maxDepth, depth, iter_pre, iter_post, iter_final)
    """
    This function is a recursive V_Cycle Multigrid for solving the Diffusion Equation on a uniform grid
    
    phi:                The matrix x in Ax=b, is either concentration (phi) or error (e)
    diffRegion:         A BitMatrix indicating per gridcell whether that cell is in the diffusion region
    bulk_value:         The boundary value of the compound outside of the diffusion region
                        so, the value in the bulk/dirichlet boundary condition
    f:                  The matrix b in Ax = b, is either rhs of diffusion-reaction equation, or residuals after residuals restriction
    L_0:                The basis of the laplacian stencil
    L_restriction:      The stencil for restriction
    L_prolongation:     The stencil for prolongation
    maxDepth:           The minimum number of gridcells in the highest recursion of the V-cycle, 
                        i.e, at the deepest level of recursion, the grid is at least (maxDepth, maxDepth)
    depth:              The current depth (starting at 0)
    iter_pre:           The number of smoothing steps before recursion
    iter_post:          Thenumber of smoothing steps after error correction
    iter_final:         The number of smoothing steps on deepest level

    returns:
    phi:                The matrix x advanced towards the solution in Ax=b
    """

    include(string(pwd(), "\\lib\\diffusion\\smoothing.jl"))
    include(string(pwd(), "\\lib\\diffusion\\prolongation.jl"))
    include(string(pwd(), "\\lib\\diffusion\\restriction.jl"))
    include(string(pwd(), "\\lib\\diffusion\\residual.jl"))
    
    # Create a correct left-hand-side stencil according to A_lhs = I - 1/2 ^2d * L_0
    L_lhs = [0 0 0; 0 1 0; 0 0 0] .- ((1/2)^(2*depth)) * L_0

    # Pre-recursion smoothing
    for itr in 1:iter_pre
        phi = smoothing(phi, f, L_lhs)                          # Smooth with the matrix A
        phi = phi .*diffRegion .+ bulk_value * .!diffRegion     # But not for the bulk region, that should be Dirichlet (or in case of error, 0)
    end

    # Compute Residual errors
    r = residual(phi, f, L_lhs)
    r = diffRegion .* r                                         # bulk_region has by definition RES=0

    # Restriction
    rhs = restriction(r, L_restriction)                         # Restrict residuals to a coarser grid. This new grid is the rhs in Ae = r
    next_diffRegion = diffRegion[1:2:end, 1:2:end]              # Adjust diffusion region similarly as the restriction has been done

    eps = zeros(size(rhs))                                      # Initial guess is 0's as when all errors are 0 we have perfect solution

    # Stop recursion at smallest grid size, otherwise continue recursion
    next_bulk_value = 0                                                         # Next iterations deal with errors and Residuals. In bulk RES and error are 0.
    if ceil(sqrt(length(phi))) <= maxDepth                                      # If the length of a dimension is the smallest size or smaller

        # stop recursion, but do some final smoothing with the last restricted residuals and error
        L_lhs_deeper = [0 0 0; 0 1 0; 0 0 0] .- ((1/2)^(2*(depth+1))) * L_0     # Create a new stencil for this final depth
        for itr in 1:iter_final
            eps = smoothing(eps, rhs, L_lhs_deeper)                             # Smooth the errors
            eps = eps .* next_diffRegion                                        # But only for diffusion region
        end
    else
        # recursively go a level deeper
        eps = V_Cycle!(eps, next_diffRegion, next_bulk_value, rhs, L_0, L_restriction, L_prolongation, maxDepth, depth+1, iter_pre, iter_post, iter_final)
    end

    # Prolongation or eps and Correction of phi
    phi = phi + prolongation(eps, L_prolongation, size(phi))

    # Post recursion smoothing
    for itr in 1:iter_post
        phi = smoothing(phi, f, L_lhs0)                         # Smooth with the matrix A
        phi = phi .* diffRegion .+ bulk_value * .!diffRegion    # But not for bulk region
    end
    return phi
end
