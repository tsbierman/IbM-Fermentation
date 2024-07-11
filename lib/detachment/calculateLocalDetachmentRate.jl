function calculateLocalDetachmentRate(i, j, kDet, grid_float, xcentre, ycentre)
    """
    This function calculates the local detachment speed with the formula Fdet = kdet * d^2

    Arguments
    i, j:               Gridcell indices along the x and y axis, respectively
    kDet:               Detachment constant determining how fast detachment is
    grid_float:         A "Float" struct containing grid parameters of type Float64
    xcentre, ycentre:   x and y coordinates of the aggregate centre

    Returns
    Fdet:               Speed of detachment
    """

    # Get centre of grid cell
    x = (i - 0.5) * grid_float.dx
    y = (j - 0.5) * grid_float.dx

    # Calculate distance from centre
    d = sqrt((x - xcentre)^2 + (y - ycentre)^2)

    # Calculate Fdet
    Fdet = kDet * d^2 # [m/h]

    # If it is the centre, correct to infinity
    if d == 0
        Fdet = Inf
    end

    return Fdet
end
