# Reproduce every figure and the summary table of the BAG paper (../draft).
# Run with threads, e.g.  julia -t auto --project=. scripts/make_paper.jl
# Basins are cached in data/paper_basins, so re-runs only redo the analysis.
# Each script runs in its own module; packages are compiled only once.
for s in ("fig_sierpinski.jl", "fig_newton.jl", "fig_duffing.jl",
          "fig_magnetic.jl", "fig_feudel.jl", "table_summary.jl")
    @info "Running $s"
    path = joinpath(@__DIR__, s)
    Core.eval(Main, :(module $(Symbol(splitext(s)[1])) include($path) end))
end
