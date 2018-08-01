# Run `julia make.jl` in this folder to generate .html pages in a build/
# directory, then open index.md for docs

push!(LOAD_PATH,"../src/")
using Documenter, DPMeansClustering

makedocs(
    format = :html,
    sitename = "DPMeansClustering.jl",
    modules = [DPMeansClustering],
    pages = [
        "index.md",
        "dp.md"
    ]
)

