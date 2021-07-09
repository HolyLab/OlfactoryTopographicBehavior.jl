module OlfactoryTopographicBehavior

using DataFrames
using Unitful: Hz, s

export session_frame, trialranges

struct TrialRanges
    trialrange::UnitRange{Int}
    soundrange::UnitRange{Int}
    lickranges::Vector{UnitRange{Int}}
end
TrialRanges(tr::AbstractUnitRange, sr::AbstractUnitRange) = TrialRanges(tr, sr, UnitRange{Int}[])

function session_frame(ai; odortiming=1, vaccuum=2, soundlick=3, camera=4, sniff=5, lickometer=6, odordirection=7, lick=8)
    return DataFrame("odortiming" => view(ai, :, odortiming),
                     "soundlick" => view(ai, :, soundlick),
                     "sniff" => view(ai, :, sniff),
                     "odordirection" => view(ai, :, odordirection),
                     "lick" => view(ai, :, lick)
                     )
end

function trialranges(df::DataFrame; fs=1000Hz, soundduration=1s, trialoffset=-1s, trialduration=10s, lohi_sound=(1, 3), lohi_lick=(0.25, 0.75), tprecision=0.1)
    soundlick, lick = df.soundlick, df.lick
    alltrs = TrialRanges[]
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
            sr = i:j-1
            tr = clamp(i + round(Int, float(trialoffset*fs)), 1, n) : clamp(i + round(Int, (trialduration + trialoffset) * float(fs)), 1, n)
            trs = TrialRanges(tr, sr)
            # Detect all licks in the trial
            j = advance_to_lo(lick, tr, lohi_lick)
            while j < last(tr)
                j = advance_to_hi(lick, j, last(tr), lohi_lick)
                j >= last(tr) && break
                k = advance_to_lo(lick, j, last(tr), lohi_lick)
                push!(trs.lickranges, j:k-1)
                j = k
            end
            push!(alltrs, trs)
            i = last(tr)
        else
            error("soundlick pulse of duration $thi starting at $i not recognized")
        end
    end
    return alltrs
end

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
