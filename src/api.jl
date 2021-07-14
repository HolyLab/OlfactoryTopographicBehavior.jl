"""
    dfr = session_frame(ai; odortiming=1, vaccuum=2, soundlick=3, camera=4, sniff=5, lickometer=6, odordirection=7, lickbinary=8, fs=$default_fs)

Construct a trial session DataFrame from an input matrix `ai`. The keyword arguments provide the column index for each signal.
The return value is a `DataFrameRate` which also encodes the sampling frequency (as supplied by `fs`).
"""
function session_frame(ai; odortiming=1, vaccuum=2, soundlick=3, camera=4, sniff=5, lickometer=6, odordirection=7, lickbinary=8, fs=default_fs)
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
    return DataFrameRate(df, fs)
end

"""
    trialdatas(dfr::DataFrameRate; fcarrierlick=60Hz, soundduration=1s, trialoffset=-1s, trialduration=10s, minlickduration=0.005s, maxlickamp=0.1, lohi_sound=(1, 3), lohi_odor=(0.5, 3.0), tprecision=0.1)
    trialdatas(df::DataFrame; fs=$default_fs, fcarrierlick=60Hz, soundduration=1s, trialoffset=-1s, trialduration=10s, minlickduration=0.005s, maxlickamp=0.1, lohi_sound=(1, 3), lohi_odor=(0.5, 3.0), tprecision=0.1)

Collect data on all trials in `df`, which should be a DataFrame with colums "soundlick", "lickometer", "odortiming", and "odordirection".
Each should be an analog signal with timing data.
Supplying `dfr` provides the DataFrame and the sampling rate.

The keyword arguments allow you to set:
- `fs`: the sampling frequency (affects conversion to physical units)
- `fcarrierlick`: the frequency of the carrier signal for lick detection (grounding this signal indicates a lick)
- `soundduration`: the expected duration of the sound pulse. This is used to find trials.
- `tprecision`: fractional "slop" in timing permitted for identifying a sound pulse.
- `trialoffset`: start of the trial relative to the sound pulse onset
- `trialduration`: length of the trial
- `minlickduration`: the minimum length of grounding the lick carrier signal to count as a lick
- `maxlickamp`: maximum allowed deviation of carrier voltage during a lick (as a fraction of lick carrier stddev)
- `lohi_sound`: voltage thresholds for sound lo and hi (a sound pulse starts when this goes hi)
- `lohi_odor`: voltage threshold for odor timing (odor pulse starts when this goes hi)
"""
function trialdatas(df::DataFrame; fs=default_fs, fcarrierlick=60Hz, soundduration=1s, trialoffset=-1s, trialduration=10s, minlickduration=0.005s, maxlickamp=0.1, lohi_sound=(1, 3), lohi_odor=(0.5, 3.0), tprecision=0.1, kwargs...)
    soundlick, lick, odortiming, dir = df.soundlick, df.lickometer, df.odortiming, df.odordirection
    checkcarrier(lick, fs, fcarrierlick; kwargs...)
    lickthresh = std(lick) * maxlickamp
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
            icue = advance_to_hi(soundlick, j, n, lohi_sound)
            iodor = advance_to_hi(odortiming, tr, lohi_odor)
            jodor = advance_to_lo(odortiming, iodor, last(tr), lohi_odor)
            oi = (iodor-first(tr))/fs .. (jodor-first(tr))/fs
            trs = TrialData(tr, si, (icue - first(tr))/fs, oi, dir[i])
            addlicks!(trs.lickonsets, view(lick, tr); fs, minlickduration, lickthresh)

            push!(alltrs, trs)
            i = last(tr)
        else
            error("soundlick pulse of duration $(uconvert(s, thi)) starting at $i not recognized")
        end
    end
    return alltrs
end
trialdatas(dfr::DataFrameRate; kwargs...) = trialdatas(dfr.df; fs=dfr.fs, kwargs...)

# Detect all licks in the trial. `lick` should be just the voltages recorded during a single trial.
function addlicks!(out, lick; fs, minlickduration, lickthresh)
    k = k1 = firstindex(lick)
    while k < lastindex(lick)
        k += 1
        abs(lick[k]) < lickthresh || continue        # if we're not close to ground, don't look for a lick
        kend = k
        while kend < lastindex(lick)
            kend += 1
            abs(lick[kend]) < lickthresh || break    # keep looking forward until an above-threshold deviation from ground
        end
        if kend - k > fs*minlickduration
            push!(out, (k - k1)/fs)                  # we stayed grounded for long enough, count it as a lick
        end
        k = kend
    end
    return out
end

"""
    is_lick(td::TrialData; deadline=4s)

Return `true` if the animal licked. Only consider licks occurring after the cue time and within `deadline` thereafter.
"""
function is_lick(td::TrialData; deadline=4s)
    tstart = td.cuetime # minimum(td.odorinterval)
    firstidx = searchsortedfirst(td.lickonsets, tstart)
    firstidx > lastindex(td.lickonsets) && return false
    return td.lickonsets[firstidx] - tstart < deadline
end
