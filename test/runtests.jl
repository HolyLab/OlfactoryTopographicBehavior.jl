using OlfactoryTopographicBehavior
using DataFrames
using Unitful: Hz, s
using InvertedIndices
using Test

@testset "OlfactoryTopographicBehavior.jl" begin
    df = DataFrame("soundlick" => [0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                   "lick"      => [0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0])
    alltrs = trialranges(df; lohi_sound = (0.25, 0.75), fs=2Hz, soundduration=2s)
    trs = only(alltrs)
    @test all(df.soundlick[trs.soundrange] .== 1)
    @test all(df.soundlick[Not(trs.soundrange)] .== 0)
    @test trs.trialrange == 2:length(df.soundlick)
    @test length(trs.lickranges) == 3
    lick = copy(df.lick)
    foreach(trs.lickranges) do lr
        @test all(lick[lr] .== 1)
        lick[lr] .= 0
    end
    @test all(iszero, lick)
    @test_throws ErrorException trialranges(df; lohi_sound = (0.25, 0.75), fs=2Hz, soundduration=3s)
end
