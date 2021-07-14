export plottrial, trialaxes

iv2t_s(iv) = (minimum(iv)/(1s), maximum(iv)/(1s))

"""
    plottrial(dfr::DataFrameRate, td::TrialData; axsniff, axtiming)
    plottrial(signal::AbstractVector, td::TrialData; fs, axsniff, axtiming)

Plot a trial in the specified axes. See [`trialaxes`](@ref) to create them.
"""
function plottrial(signal::AbstractVector, td::TrialData; fs, axsniff, axtiming, showlegend=true)
    axsniff.cla()
    axtiming.cla()
    sniff = signal[td.trialrange]
    t = (0:length(sniff)-1)/fs
    axsniff.plot(t * (1Hz), sniff, color="green")
    axsniff.set_xlabel("Time (s)")
    axsniff.set_ylabel("Chest expansion (V)")
    xl = axsniff.get_xlim()
    axtiming.plot(iv2t_s(td.soundinterval), (1, 1), color="blue")
    tcue = convert(Float64, td.cuetime / (1s))
    axtiming.plot(tcue, 1.5, "kv")
    axtiming.plot(iv2t_s(td.odorinterval), (2, 2), color="cyan")
    axtiming.text(iv2t_s(td.odorinterval)[end], 2, string(td.odordirection) * (is_lick(td) ? "+" : "-"))
    tlick = td.lickonsets / (1s)
    if !isempty(tlick)
        axtiming.plot([tlick'; tlick'], [fill(2.75, length(tlick))'; fill(3.25, length(tlick))'], color="red")
    end
    axtiming.set_xlim(xl)
    axtiming.set_ylim((0.5, 3.75))
    axtiming.set_xticks([])
    axtiming.set_yticks([])
    axsniff.set_frame_on(false)
    axtiming.set_frame_on(false)
    axtiming.invert_yaxis()
    if showlegend
        axtiming.legend(("sound", "cue", "odor", "licks"))
    end
    return nothing
end

plottrial(dfr::DataFrameRate, td::TrialData; kwargs...) = plottrial(dfr.sniff, td; fs=dfr.fs, kwargs...)


"""
    axsniff, axtiming = plottrial(dfr::DataFrameRate, td::TrialData, fig::PyPlot.Figure)

Plot a trial in the given `fig`ure.
"""
function plottrial(dfr::DataFrameRate, td::TrialData, fig::PyPlot.Figure)
    axsniff, axtiming = trialaxes(fig)
    plottrial(dfr, td; axsniff, axtiming)
    return axsniff, axtiming
end

"""
    axsniff, axtiming = trialaxes(fig)

Create axes for plotting trials.
"""
function trialaxes(fig::PyPlot.Figure)
    gs = PyPlot.plt.GridSpec(2, 1, fig, height_ratios=[1, 3])
    return fig.add_subplot(gs[2]), fig.add_subplot(gs[1])
end
