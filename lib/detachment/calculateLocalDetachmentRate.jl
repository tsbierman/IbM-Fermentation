function calculateLocalDetachmentRate(i, j, kDet, grid, xcentre, ycentre)
    """
    Calculate the local detachment speed: Fdet = kdet * d^2
    i, j: gridcell indices along the x and y axis respectively
    kDet: Detachment constant determining how fast detachment is
    grid: struct containing all information regarding spacial discretization
    xcentre, ycentre, x and y centre of the coordinates

    Returns:
    Fdet: speed of detachment
    """

    # Get centre of grid cell
    x = (i - 0.5) * grid.dx
    y = (j - 0.5) * grid.dx

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
