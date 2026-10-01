module Plots

import ..Data

function plotValsOnMap end
function plotValsOnMap! end
function plotMapGrid! end
function plotZonalMean end
function plotZonalMean! end
function plotAMOC end
function plotEnsembleSpread end
function plotTimeseries end
function plotTimeseries! end
function plotTempGraph end
function makeScatterPlot end
function plotECDF end
function plotPDF end
function plotExpectedECS end


function savePlot end
function plotDistMatrices end
function convertKgsToSv! end
function mapValsToColors end
function addMinorGrid! end
# function makeSubplots end 

function plotWeights end
function plotDistances end
function plotDistancesIndependence end
function plotCRPSSPseudoObs end
function crpssBoxPlot end
function boxplotMCMCWeights end
function densityMCMCWeights end
function plotCorrWeights end


# All functions above are only stubs, their methods are defined in the package
# extension MakieExt, which is only loaded when CairoMakie and GeoMakie are loaded.
# Add a hint to MethodErrors thrown by these stubs if the extension isn't loaded.
function __init__()
    Base.Experimental.register_error_hint(MethodError) do io, exc, argtypes, kwargs
        f = exc.f
        if f isa Function && parentmodule(f) === Plots &&
           isnothing(Base.get_extension(parentmodule(Plots), :MakieExt))
            printstyled(io, "\n\nModelWeights.Plots.$(nameof(f)) is only available when " *
                "the plotting extension is loaded. Run `using CairoMakie, GeoMakie` first.";
                color = :cyan)
        end
    end
end

end