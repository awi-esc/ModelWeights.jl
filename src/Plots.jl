module Plots

import StatsBase.ecdf

using CairoMakie
using ColorSchemes
using Dates
using DimensionalData
using GeoMakie
using Statistics
using TextWrap
using YAXArrays

CairoMakie.activate!(type = "svg")


include("plot-utils.jl")
include("plot-data.jl")
include("plot-weights.jl")

using ..Data

end