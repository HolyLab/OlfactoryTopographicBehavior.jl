const TU = typeof(1.0f0s)
const TI = typeof(1.0f0s .. 2.0f0s)
const default_fs = 1000.0f0Hz

struct TrialData
    trialrange::UnitRange{Int}
    soundinterval::TI
    cuetime::TU
    odorinterval::TI
    odordirection::Int
    lickonsets::Vector{TU}
end
TrialData(tr::AbstractUnitRange, si::AbstractInterval, ct::Quantity, oi::AbstractInterval, od) = TrialData(tr, si, ct, oi, od, TU[])

function Base.show(io::IO, td::TrialData)
    println(io, "Trial over indices ", td.trialrange, ':')
    println(io, "  sound: ", td.soundinterval)
    println(io, "  cue: ", td.cuetime)
    println(io, "  odor delivery: ", td.odorinterval)
    println(io, "  odor direction: ", td.odordirection)
    print(io, "  lick times: ")
    Base.show_delim_array(io, td.lickonsets, '[', ",", ']', false)  # skips printing of unit
end

struct DataFrameRate
    df::DataFrame
    fs::typeof(default_fs)
end

function Base.show(io::IO, dfr::DataFrameRate)
    show(io, dfr.df)
    print(io, "\nRecorded at ", dfr.fs)
end

function Base.getproperty(dfr::DataFrameRate, name::Symbol)
    name === :fs && return getfield(dfr, :fs)
    df = getfield(dfr, :df)
    name === :df && return df
    return getproperty(df, name)
end
