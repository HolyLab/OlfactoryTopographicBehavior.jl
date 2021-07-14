using OlfactoryTopographicBehavior
using DataFrames
using Unitful: Hz, s
using IntervalSets
using Test

@testset "OlfactoryTopographicBehavior.jl" begin
    # TrialData extraction
    df = DataFrame("soundlick"     => [0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                   "lickometer"    => [1,-1, 1,-1, 1, 0, 0, 0, 1,-1, 1, 0, 1,-1, 1],
                   "odortiming"    => [0, 0, 0, 0, 0, 0, 0, 0, 4, 4, 0, 0, 0, 0, 0],
                   "odordirection" => [1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1])
    fs = 2Hz
    tds = trialdatas(df; lohi_sound = (0.25, 0.75), fs, soundduration=2s, fcarrierlick=1Hz, minlickduration=1s)
    td = only(tds)
    @test td.trialrange == 2:length(df.soundlick)
    @test td.soundinterval == 1s .. 3s
    @test only(td.lickonsets) == 2s
    @test td.odordirection == 2
    @test td.odorinterval == 3.5s .. 4.5s
    tds = trialdatas(df; lohi_sound = (0.25, 0.75), fs, soundduration=2s, fcarrierlick=1Hz, minlickduration=0s)
    td = only(tds)
    @test td.lickonsets == [2s, 5s]
    @test_throws ErrorException("peak frequency 1.0 Hz does not match expected 60 Hz") trialdatas(df; lohi_sound = (0.25, 0.75), fs, soundduration=3s)
    @test_throws ErrorException("soundlick pulse of duration 2.0 s starting at 4 not recognized") trialdatas(df; lohi_sound = (0.25, 0.75), fs, soundduration=3s, fcarrierlick=1Hz, minlickduration=1s)

    # Display
    io = IOBuffer()
    print(io, td)
    str = String(take!(io))
    @test startswith(str, "Trial over indices 2:15")

    # Construction from a matrix
    mat = [df.soundlick df.lickometer df.odortiming df.odordirection]
    dfmat = session_frame(mat; fs, soundlick=1, lickometer=2, odortiming=3, odordirection=4)
    @test dfmat.df == df
    print(io, dfmat)
    str = String(take!(io))
    @test endswith(str, "Recorded at 2.0f0 Hz")
    @test dfmat.lickometer == df.lickometer
    @test trialdatas(dfmat;  lohi_sound = (0.25, 0.75), soundduration=2s, fcarrierlick=1Hz, minlickduration=1s) ==
          trialdatas(df; fs, lohi_sound = (0.25, 0.75), soundduration=2s, fcarrierlick=1Hz, minlickduration=1s)

    # lick detection
    td = OlfactoryTopographicBehavior.TrialData(1:0, 1s .. 2s, 1s, 1.5s .. 3.0s, -1, [0.8s])
    @test !is_lick(td)
    push!(td.lickonsets, 6s)
    @test !is_lick(td)
    @test  is_lick(td; deadline=10s)
    empty!(td.lickonsets)
    @test !is_lick(td)
    push!(td.lickonsets, 1.2s)
    @test  is_lick(td)
    empty!(td.lickonsets)
    push!(td.lickonsets, 0.8s, 1.2s)
    @test  is_lick(td)
end
