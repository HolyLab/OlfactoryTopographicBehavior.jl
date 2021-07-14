module OlfactoryTopographicBehavior

using IntervalSets: minimum
using DataFrames
using IntervalSets
using Unitful: Hz, s, uconvert, Quantity
using DSP
using Statistics
using Requires

export TrialData, collect_by_trialtype, is_lick, session_frame, trialdatas

include("types.jl")
include("utils.jl")
include("api.jl")

## __init__

function __init__()
    @require PyPlot = "d330b81b-6aea-500a-939a-2ce795aea3ee" include("plotting.jl")
end

end
