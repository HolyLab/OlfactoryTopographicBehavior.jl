sniffslope(df::DataFrameRate, σ::Quantity) = slope(smooth(df.sniff, σ; fs=df.fs))

function sniffgrams(sniff, tds::AbstractVector{TrialData}; fs=default_fs, sniffwindow=0s..2s, Δf=1/width(sniffwindow))
    sgs = Matrix{eltype(sniff)}[]
    freq = time = nothing
    n = ceil(Int, convert(AbstractFloat, fs/Δf))
    for td in tds
        ostart = minimum(td.odorinterval)
        tr = first(td.trialrange) .+ iv2r(ostart + minimum(sniffwindow) .. ostart + maximum(sniffwindow), fs)
        if length(tr) < n
            tr = first(tr) : first(tr) + n -1
        end
        if last(tr) <= lastindex(sniff)
            sg = Periodograms.spectrogram(view(sniff, tr), n; fs=fs/(1Hz))
            freq, time = sg.freq, sg.time
            push!(sgs, sg.power)
        end
    end
    nt = length(time)
    Δt = width(sniffwindow)/nt
    return range(minimum(sniffwindow)+Δt/2, stop=maximum(sniffwindow)-Δt/2, length=nt), freq*Hz, sgs
end
sniffgrams(dfr::DataFrameRate, tds::AbstractVector{TrialData}; kwargs...) = sniffgrams(dfr.df.sniff, tds; fs=dfr.fs, kwargs...)

function sniffgrams_power(freq, sgs; finterval=4Hz .. 10Hz)
    idx = searchsortedfirst(freq, minimum(finterval)) : searchsortedfirst(freq, maximum(finterval))
    return [sum(sg[idx,:]) for sg in sgs]
end
