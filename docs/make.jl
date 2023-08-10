using EulerLagrange
using Documenter

makedocs(;
    modules=[EulerLagrange],
    authors="Michael Kraus",
    repo="https://github.com/JuliaGNI/EulerLagrange.jl/blob/{commit}{path}#L{line}",
    sitename="EulerLagrange.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JuliaGNI.github.io/EulerLagrange.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo   = "github.com/JuliaGNI/EulerLagrange.jl",
    devurl = "latest",
    devbranch = "main",
)
