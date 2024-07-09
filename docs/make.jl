using Documenter, PANDA

makedocs(
    sitename="PANDA.jl",
    modules=[PANDA],
    pages=[
        "Home" => "index.md",
        "ClaDS" => ["Manual" => "clads/tutorial.md"]
    ],
    authors="Helene Morlon, Odile Maliet",
    checkdocs = :none,
)

deploydocs(
    repo="github.com/hmorlon/PANDA.jl",
)
