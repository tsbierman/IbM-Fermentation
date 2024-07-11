function update_bacterial_radius!(bac_vecfloat, constants_float)
    """
    This function updates all bacterial radii based on their masses.

    Arguments
    bac_vecfloat:           A "VectorFloat" struct containing bacterial parameters of type Vector{Float64}
    constants_float:        A "Float" struct containing simulation constants of type Float64

    Returns
    bac_vecfloat:           A "VectorFloat" struct containing bacterial parameters of type Vector{Float64} with updated radius
    """
    
    # Determine radius based on the mass and spherical cell
    bac_vecfloat.radius = ((bac_vecfloat.molarMass .* constants_float.bac_MW ./ constants_float.bac_rho) .* (3 ./ (4 .* pi))) .^ (1 ./ 3)
    return bac_vecfloat
end
