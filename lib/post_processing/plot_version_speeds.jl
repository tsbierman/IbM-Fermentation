function HEX2RGB(c)

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
    "Makes a plot of the computation speeds of the different versions"

    file = XLSX.readxlsx(string(pwd(), "\\results\\Speed Results.xlsx"))
    data = file["Sheet2"]["A6:I11"]
    sim_numbers = data[:,1]
    versions = file["Sheet2"]["B6:B11"]

    shortened_data = convert(Matrix{Float32}, data[:,3:end])

    plotting_data = shortened_data[:,1:5]
    v102_outlier = shortened_data[2, 6]
    v200_outliers = shortened_data[3, 6:7]
    v201_outlier = shortened_data[4, 6]
    v202_outliers = shortened_data[5,6:7]

    means = mean(plotting_data, dims=2)
    standard_deviations = std(plotting_data, dims=2)

    # rC, gC, bC = HEX2RGB(["#88A3EE", "#163CA5"])
    rC, gC, bC = HEX2RGB(["#B8E186", "#276419"])
    my_blue = RGB(rC[1]/255, gC[1]/255, bC[1]/255)
    my_dark_blue = RGB(rC[2]/255, gC[2]/255, bC[2]/255)

    bar(means,
        xticks=(1:length(means), versions),
        xlabel= "Version of IbM",
        ylabel= "Time (s)",
        c = my_blue,
        ylim = [0,2000],
        yerror = standard_deviations,
        legend = false,
        dpi=300)

    scatter!([2], [v102_outlier], c= my_dark_blue)
    scatter!([3,3], v200_outliers, c= my_dark_blue)
    scatter!([4], [v201_outlier], c= my_dark_blue)
    scatter!([5,5], v202_outliers, c= my_dark_blue)

    savefig(string(pwd(), "\\results\\Version_speeds.png"))

end
