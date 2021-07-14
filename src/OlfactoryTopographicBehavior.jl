module OlfactoryTopographicBehavior

using IntervalSets: minimum
using DataFrames
using IntervalSets
using Unitful: Hz, s, uconvert, Quantity
using DSP
using Statistics

export TrialData, is_lick, session_frame, trialdatas

include("types.jl")
include("utils.jl")
include("api.jl")

end
