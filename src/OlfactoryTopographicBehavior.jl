module OlfactoryTopographicBehavior

using DataFrames
using IntervalSets
using Unitful: Hz, s, uconvert

export is_lick, session_frame, trialdatas

const TU = typeof(1.0s)
const TI = typeof(1.0s .. 2.0s)

struct TrialData
    trialrange::UnitRange{Int}
    soundinterval::TI
    odorinterval::TI
    odordirection::Int
    lickonsets::Vector{TU}
end
TrialData(tr::AbstractUnitRange, si::AbstractInterval, oi::AbstractInterval, od) = TrialData(tr, si, oi, od, TU[])

function Base.show(io::IO, td::TrialData)
    println(io, "Trial over indices ", td.trialrange, ':')
    println(io, "  sound: ", td.soundinterval)
    println(io, "  odor delivery: ", td.odorinterval)
    println(io, "  odor direction: ", td.odordirection)
    print(io, "  lick times: ")
    Base.show_delim_array(io, td.lickonsets, '[', ",", ']', false)  # skips printing of unit
end

"""
    df = session_frame(ai; odortiming=1, vaccuum=2, soundlick=3, camera=4, sniff=5, lickometer=6, odordirection=7, lick=8)

Construct a trial session DataFrame from an input matrix `ai`. The keyword arguments provide the column index for each signal.
"""
function session_frame(ai; odortiming=1, vaccuum=2, soundlick=3, camera=4, sniff=5, lickometer=6, odordirection=7, lick=8)
    df = DataFrame()
    if soundlick ∈ axes(ai, 2)
        df.soundlick = view(ai, :, soundlick)
    end
    if lickometer ∈ axes(ai, 2)
        df.lickometer = view(ai, :, lickometer)
    end
    if odortiming ∈ axes(ai, 2)
        df.odortiming = view(ai, :, odortiming)
    end
    if odordirection ∈ axes(ai, 2)
        df.odordirection = view(ai, :, odordirection)
    end
    if sniff ∈ axes(ai, 2)
        df.sniff = view(ai, :, sniff)
    end
    return df
end

"""
    trialdatas(df::DataFrame; fs=1000Hz, soundduration=1s, trialoffset=-1s, trialduration=10s, lohi_sound=(1, 3), lohi_lick=(-0.2, -0.05), lohi_odor=(0.5, 3.0), tprecision=0.1)

Collect data on all trials in `df`, which should be a DataFrame with colums "soundlick", "lickometer", "odortiming", and "odordirection".
Each should be an analog signal with timing data.

The keyword arguments allow you to set:
- `fs`: the sampling frequency (affects conversion to physical units)
- `soundduration`: the expected duration of the sound pulse. This is used to find trials
- `tprecision`: fractional "slop" in timing permitted for identifying a sound pulse.
- `trialoffset`: start of the trial relative to the sound pulse onset
- `trialduration`: length of the trial
- `lohi_sound`: voltage thresholds for sound lo and hi (a sound pulse starts when this goes hi)
- `lohi_lick`: voltage thresholds for lick lo and hi (negative deflections correspond to licks)
- `lohi_odor`: voltage threshold for odor timing (odor pulse starts when this goes hi)
"""
function trialdatas(df::DataFrame; fs=1000Hz, soundduration=1s, trialoffset=-1s, trialduration=10s, lohi_sound=(1, 3), lohi_lick=(-0.2, -0.05), lohi_odor=(0.5, 3.0), tprecision=0.1)
    soundlick, lick, odortiming, dir = df.soundlick, df.lickometer, df.odortiming, df.odordirection
    alltrs = TrialData[]
    i, n = 1, length(soundlick)
    i = advance_to_lo(soundlick, i, n, lohi_sound)
    while i < n
        i = advance_to_hi(soundlick, i, n, lohi_sound)
        i >= n && break
        # We're at a rising edge. When does it fall again?
        j = advance_to_lo(soundlick, i, n, lohi_sound)
        thi = (j-i)/fs   # duration, in time units, of the hi pulse
        if (1-tprecision)*soundduration <= thi <= (1+tprecision)*soundduration
            # We're in a sound marker. Initiate a trial.
            tr = clamp(i + round(Int, float(trialoffset*fs)), 1, n) : clamp(i + round(Int, (trialduration + trialoffset) * float(fs)), 1, n)
            si = (i-first(tr))/fs .. (j-first(tr))/fs
            iodor = advance_to_hi(odortiming, tr, lohi_odor)
            jodor = advance_to_lo(odortiming, iodor, last(tr), lohi_odor)
            oi = (iodor-first(tr))/fs .. (jodor-first(tr))/fs
            trs = TrialData(tr, si, oi, dir[i])
            # Detect all licks in the trial
            j = advance_to_hi(lick, tr, lohi_lick)   # lick signal is hi when not licking
            while j < last(tr)
                j = advance_to_lo(lick, j, last(tr), lohi_lick)
                j >= last(tr) && break
                push!(trs.lickonsets, (j - first(tr))/fs)
                j = advance_to_hi(lick, j, last(tr), lohi_lick)
            end
            push!(alltrs, trs)
            i = last(tr)
        else
            error("soundlick pulse of duration $(uconvert(s, thi)) starting at $i not recognized")
        end
    end
    return alltrs
end

"""
    is_lick(td::TrialData; deadline=4s)

Return `true` if the animal licked. Only consider licks occurring after the beginning of the sound pulse and within `deadline` thereafter.
"""
function is_lick(td::TrialData; deadline=4s)
    tstart = minimum(td.soundinterval)
    firstidx = searchsortedfirst(td.lickonsets, tstart)
    firstidx > lastindex(td.lickonsets) && return false
    return td.lickonsets[firstidx] - tstart < deadline
end

## Utilities

islo(trace, idx, lohi) = trace[idx] < lohi[1]
ishi(trace, idx, lohi) = trace[idx] > lohi[2]

function advance_to(isval::Function, trace, idxstart::Integer, idxstop::Integer, lohi)
    i = idxstart
    while i < idxstop && !isval(trace, i, lohi)
        i += 1
    end
    return i
end
advance_to(isval::Function, trace, startstop::AbstractUnitRange{Int}, lohi) =
    advance_to(isval, trace, first(startstop), last(startstop), lohi)

advance_to_lo(args...) = advance_to(islo, args...)
advance_to_hi(args...) = advance_to(ishi, args...)

end
