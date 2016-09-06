module PyPlotHelper

export
showmatrix,
manhattanplot,
zhist,
probhist,
uniform_quantile,
zz_scratterplot

using PyCall
using PyPlot
using Distributions

include("plots.jl")

end
