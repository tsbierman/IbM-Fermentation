function HEX2RGB(c)
    """
    Function to convert hexcode to RGB
    """
    rC, gC, bC = zeros(length(c)), zeros(length(c)), zeros(length(c))
    for idx in eachindex(c)
        color = c[idx]
        RGB = parse(UInt, color[2:end], base=16)
        rC[idx] = (RGB >> 16) & 0xFF
        gC[idx] = (RGB >> 8) & 0xFF
        bC[idx] = RGB & 0xFF
    end
    return rC, gC, bC
end

function plot_modelparts_speeds()
    """
    Function to plots the time the model spends in the several parts of the code
    """

    # Select and get data
    file = XLSX.readxlsx(string(pwd(), "\\results\\Model Parts Times.xlsx"))
    names = file["Julia"]["A10:A20"]
    julia_data = file["Julia"]["C2:M8"]
    julia2_data = file["Julia2"]["C2:M5"]
    matlab_data = file["MATLAB"]["C2:M5"]

    # Get names and convert data
    names = convert(Matrix{String}, names)
    julia_data = convert(Matrix{Float32}, julia_data)
    julia2_data = convert(Matrix{Float32}, julia2_data)
    matlab_data = convert(Matrix{Float32}, matlab_data)

    # Colors used
    c = ["#276419", "#4D9221", "#7FBC41", "#B8E186", "#E6F5D0", "#F7F7F7", "#FDE0EF", "#8E0152", "#F1B6DA", "#C51B7D", "#DE77AE"]
    rC, gC, bC = HEX2RGB(c)
    colors_database = [RGB(rC[i]/255, gC[i]/255, bC[i]/255) for i in eachindex(rC)]

    # Calculate average
    julia_average = mean(transpose(julia_data), dims=2) ./ 100
    julia2_average = mean(transpose(julia2_data), dims=2) ./ 100
    matlab_average = mean(transpose(matlab_data), dims=2) ./ 100

    # Plot pie plots for julia and MATLAB
    pie(names, julia_average, 
        legendposition=:outerbottom,
        legend_column = 3,
        legend_font_pointsize = 7,
        legend = false,
        c = colors_database,
        title="Julia v2.0.2",
        dpi=300)

    savefig(string(pwd(), "\\results\\julia_parts_speeds.png"))

    pie(names, julia2_average, 
    legendposition=:outerbottom,
    legend_column = 3,
    legend_font_pointsize = 7,
    legend = false,
    c = colors_database,
    title="Julia v2.0.3",
    dpi=300)

    savefig(string(pwd(), "\\results\\julia2_parts_speeds.png"))

    pie(names, matlab_average, 
    legendposition=:outerbottom,
    legend_column = 3,
    legend_font_pointsize = 7,
    legend=false,
    c = colors_database,
    title="MATLAB",
    dpi=300)

    savefig(string(pwd(), "\\results\\matlab_parts_speeds.png"))
end
