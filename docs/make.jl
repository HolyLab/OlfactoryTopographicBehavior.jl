using OlfactoryTopographicBehavior
using Documenter

DocMeta.setdocmeta!(OlfactoryTopographicBehavior, :DocTestSetup, :(using OlfactoryTopographicBehavior); recursive=true)

makedocs(;
    modules=[OlfactoryTopographicBehavior],
    authors="Tim Holy <tim.holy@gmail.com> and contributors",
    repo="https://github.com/HolyLab/OlfactoryTopographicBehavior.jl/blob/{commit}{path}#{line}",
    sitename="OlfactoryTopographicBehavior.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://HolyLab.github.io/OlfactoryTopographicBehavior.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/HolyLab/OlfactoryTopographicBehavior.jl",
    devbranch="main",
)
