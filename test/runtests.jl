using OlfactoryTopographicBehavior
using DataFrames
using Unitful: Hz, s
using IntervalSets
using Test

@testset "OlfactoryTopographicBehavior.jl" begin
    # TrialData extraction
    df = DataFrame("soundlick"     => [0, 0, 0, 1, 1, 1,    1,    0, 0,     0, 0, 0, 0,     0,     0],
                   "lickometer"    => [0, 0, 0, 0, 0, -0.3, -0.1, 0, -0.25, 0, 0, 0, -0.22, -0.03, 0],
                   "odortiming"    => [0, 0, 0, 0, 0, 0,    0,    0, 4,     4, 0, 0, 0,     0,     0],
                   "odordirection" => [1, 2, 2, 2, 2, 2,    2,    2, 2,     2, 2, 1, 1,     1,     1])
    fs = 2Hz
    tds = trialdatas(df; lohi_sound = (0.25, 0.75), fs, soundduration=2s)
    td = only(tds)
    @test td.trialrange == 2:length(df.soundlick)
    @test td.soundinterval == 1s .. 3s
    @test length(td.lickonsets) == 3
    foreach(td.lickonsets) do t
        idx = Int(t*fs + first(td.trialrange))
        @test df.lickometer[idx] < -0.2
    end
    @test td.odordirection == 2
    @test td.odorinterval == 3.5s .. 4.5s
    @test_throws ErrorException("soundlick pulse of duration 2.0 s starting at 4 not recognized") trialdatas(df; lohi_sound = (0.25, 0.75), fs, soundduration=3s)

    # Display
    io = IOBuffer()
    print(io, td)
    str = String(take!(io))
    @test startswith(str, "Trial over indices 2:15")

    # Construction from a matrix
    mat = [df.soundlick df.lickometer df.odortiming df.odordirection]
    dfmat = session_frame(mat; soundlick=1, lickometer=2, odortiming=3, odordirection=4)
    @test dfmat == df
    @test all(iszero, lick)
    @test_throws ErrorException trialranges(df; lohi_sound = (0.25, 0.75), fs=2Hz, soundduration=3s)
end
