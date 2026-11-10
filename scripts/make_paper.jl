# Run with threads:  julia -t auto --project=. scripts/make_paper.jl
# Basins are cached in data/paper_basins, so re-runs only redo the analysis.
for s in ("fig_newton.jl", "paper_korcak.jl")
    @info "Running $s"
    path = joinpath(@__DIR__, s)
    Core.eval(Main, :(module $(Symbol(splitext(s)[1])) include($path) end))
end
