using StatsPlots
using Plots.PlotMeasures
filename = string(pwd(), "\\results\\3.0_O2_mass_analysis.xlsx")
file = XLSX.readxlsx(filename)

fractions = file["Masses"]["B13:G15"]
names = file["Masses"]["B12:G12"]
fractions = convert(Matrix{Float32}, fractions)

colors =["#ECE239", "#1C8C6F", "#51A9DC", "#006AAE", "#E49F07", "#D15D0A"]

groupedbar(fractions, 
           bar_position = :stack,
           bar_width=0.4,
           xticks=(1:length(fractions[:,1]), 1:length(fractions[:,1])),
           yticks=(0:0.1:1.1 , ["$(i)%" for i in 0:10:110]),
           xlabel="Simulation number",
           ylabel="Relative abundance (wt. %)",
           legend = :outerbottomright,
           linecolor = :match,
           topmargin=5mm,
           color = ["#ECE239" "#1C8C6F" "#51A9DC" "#006AAE" "#E49F07" "#D15D0A"; "#ECE239" "#1C8C6F" "#51A9DC" "#006AAE" "#E49F07" "#D15D0A"; "#ECE239" "#1C8C6F" "#51A9DC" "#006AAE" "#E49F07" "#D15D0A" ],
           label=names)

ratios = round.(fractions[:,end] ./ (1 .- fractions[:,end]), digits=0)
ratios_string = ["$(num):1" for num in ratios]

ratio1 = ratios_string[1]
ratio2 = ratios_string[2]
ratio3 = ratios_string[3]


annotate!(1, 1:3, (ratio1,10), :bottom)
annotate!(2, 1:3, (ratio2,10), :bottom)
annotate!(3, 1:3, (ratio3,10), :bottom)

savefig(string(pwd(),"\\results\\3.0_O2_mass_analysis.png"))
