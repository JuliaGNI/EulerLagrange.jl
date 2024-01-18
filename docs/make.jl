using EulerLagrange
using Documenter

makedocs(;
    modules=[EulerLagrange],
    authors="Michael Kraus",
    sitename="EulerLagrange.jl",
    format=Documenter.HTML(;
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://JuliaGNI.github.io/EulerLagrange.jl",
        assets = String[],
    ),
    pages=[
        "Home" => "index.md",
        "Hamiltonian Systems" => "hamiltonian.md",
        "Lagrangian Systems" => "lagrangian.md",
        "Degenerate Lagrangian Systems" => "degenerate_lagrangian.md",
        "Caveats" => "caveats.md",
        "Library" => "library.md",
    ],
)

deploydocs(;
    repo   = "github.com/JuliaGNI/EulerLagrange.jl",
    devurl = "latest",
    devbranch = "main",
)
