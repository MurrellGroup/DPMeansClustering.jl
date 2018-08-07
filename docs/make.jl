# Run `julia make.jl` in this folder to generate .html pages in a build/
# directory, then open index.md for docs

using Documenter, DPMeansClustering
    deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo   = "github.com/MurrellGroup/DPMeansClustering.jl.git",
    julia  = "nightly",
    osname = "osx"
makedocs()

