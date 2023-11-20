# Very important plotting of circles
using Plots
function circleShape(h, k, r)
    thet = LinRange(0, 2*pi, 500)
    return h .+ r * sin.(thet), k .+ r*cos.(thet)
end

plot(circleShape(0,0,1), seriestype = [:shape,],
legend = false, aspect_ratio=1, fillalpha = 0.2,
linealpha=0.2, c=:blue, linecolor=:blue)