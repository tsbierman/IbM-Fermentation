function granule_density!(bac_vecfloat, constants_float)
    """
    Updates the density of the granule

    Arguments
    bac_vecfloat:       A "VectorFloat" struct containing bacterial parameters of type Vector{Float64}
    constants_float:    A "Float" struct containing simulation constants of type Float64

    Returns
    bac_vecfloat:       A "VectorFloat" struct containing bacterial parameters of type Vector{Float64} with updated biological density of the granule
    """

    bac_m = bac_vecfloat.molarMass .* constants_float.bac_MW  # [g] for every bacteria

    # Approximated as a square due to irregular growth directions
    bac_vecfloat.bac_rho_bio = sum(bac_m) ./ ((maximum(bac_vecfloat.y) .- minimum(bac_vecfloat.y)) .* (maximum(bac_vecfloat.x) .- minimum(bac_vecfloat.x)) .* 2e-6)     # [g/m3] 2e-6 is the depth dimension (z-axis)
    return bac_vecfloat
end
