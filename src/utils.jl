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

function carrierfreq(lick, fs; twindow=4s, fwindow=2)
    # Check for the lickometer carrier signal
    sg = Periodograms.spectrogram(lick, nextpow(2, convert(AbstractFloat, twindow*fs)); fs=fs/(1Hz))
    pwr = vec(sum(sg.power; dims=2))[begin+1:end]   # drop the zero-frequency bin
    frq = sg.freq[begin+1:end]
    idx = argmax(pwr)
    pwrc = copy(pwr)
    pwrc[max(idx-fwindow, begin):min(idx+fwindow, end)] .= 0
    return frq[idx] * Hz, pwr[idx] / maximum(pwrc)
end

function checkcarrier(lick, fs, fcarrierlick=60Hz; frtol=0.05, powthresh=5)
    fpeak, prat = carrierfreq(lick, fs)
    abs(fpeak - fcarrierlick) < frtol * fcarrierlick || error("peak frequency $fpeak does not match expected $fcarrierlick")
    prat > powthresh || error("power ratio $prat does not meet SNR requirements")
    return nothing
end
