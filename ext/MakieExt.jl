module MakieExt

using CairoMakie, GeoMakie

import StatsBase.ecdf
import ModelWeights: Data
using ColorSchemes, Dates, DimensionalData, Distributions, Statistics, YAXArrays

import ModelWeights.Plots: plotValsOnMap, plotValsOnMap!, plotMapGrid!, plotZonalMean, plotZonalMean!, plotAMOC, plotEnsembleSpread
import ModelWeights.Plots: plotTimeseries, plotTimeseries!, plotTempGraph, makeScatterPlot, plotECDF, plotPDF, plotExpectedECS

import ModelWeights.Plots: savePlot, plotDistMatrices, convertKgsToSV!
import ModelWeights.Plots: plotWeights, plotDistances, plotDistancesIndependence, plotCRPSSPseudoObs, crpssBoxPlot, boxplotMCMCWeights, densityMCMCWeights, plotCorrWeights 

function __init__()
    CairMakie.activate!(type = "svg")
end

include("plot-utils.jl")
include("plot-data.jl")
include("plot-weights.jl")

end