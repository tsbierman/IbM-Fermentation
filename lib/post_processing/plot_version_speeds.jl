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

function plot_version_speeds()
    """
    Function to plot the speeds of the several IbM Versions"""

    # Select data
    file = XLSX.readxlsx(string(pwd(), "\\results\\Speed Results.xlsx"))
    data = file["Sheet2"]["A6:I11"]
    sim_numbers = data[:,1]
    versions = file["Sheet2"]["B6:B11"]

    # Only select required data
    shortened_data = convert(Matrix{Float32}, data[:,3:end])

    # Divide data is several parts
    MATLAB_data = shortened_data[1,1:5]
    v102_data = shortened_data[2,1:6]
    v200_data = shortened_data[3,1:7]
    v201_data = shortened_data[4,1:6]
    v202_data = shortened_data[5,1:7]
    v203_data = shortened_data[6,1:5]

    # Storage of means and standard deviations
    means = zeros(6)
    standard_deviations = zeros(6)

    # Calculate means and standard deviations
    means[1], standard_deviations[1] = mean(MATLAB_data), std(MATLAB_data)
    means[2], standard_deviations[2] = mean(v102_data), std(v102_data)
    means[3], standard_deviations[3] = mean(v200_data), std(v200_data)
    means[4], standard_deviations[4] = mean(v201_data), std(v201_data)
    means[5], standard_deviations[5] = mean(v202_data), std(v202_data)
    means[6], standard_deviations[6] = mean(v203_data), std(v203_data)

    # rC, gC, bC = HEX2RGB(["#88A3EE", "#163CA5"])
    rC, gC, bC = HEX2RGB(["#B8E186", "#276419"])
    my_blue = RGB(rC[1]/255, gC[1]/255, bC[1]/255)
    my_dark_blue = RGB(rC[2]/255, gC[2]/255, bC[2]/255)

    bar(means[1:6],
        xticks=(1:length(means), versions),
        xlabel= "Version of IbM",
        ylabel= "Time (s)",
        xlim=(0.3,6.7),
        c = my_blue,
        ylim = [0,2000],
        yerror = standard_deviations[1:6],
        legend = false,
        dpi=300)

    savefig(string(pwd(), "\\results\\Version_speeds_greens_6.png"))
end
