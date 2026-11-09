# AnalysisBasinsNets

This code base is using the [Julia Language](https://julialang.org/) and
[DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> AnalysisBasinsNets

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

You may notice that most scripts start with the commands:
```julia
using DrWatson
@quickactivate "AnalysisBasinsNets"
```
which auto-activate the project and enable local path handling from DrWatson.

## Reproducing the figures of the BAG paper

The draft lives in `../draft` (LaTeX). All its figures and Table 1 are produced by

```
julia -t auto --project=. scripts/make_paper.jl
```

which runs, in order, `scripts/fig_sierpinski.jl`, `fig_newton.jl`, `fig_duffing.jl`,
`fig_magnetic.jl`, `fig_feudel.jl` and `table_summary.jl`; each script can also be run
on its own. Figures are written to `../draft/figures` and the table to
`../draft/tables/table1.tex`.

* `src/paper_systems.jl` defines the systems with the parameters of the paper. The
  attractors are located once, then every resolution is classified with a threaded
  `AttractorsViaProximity` mapper, so labels are consistent across resolutions.
  Basins are cached in `data/paper_basins` (delete a file to recompute it).
* `src/analyze_nets.jl` builds the basin adjacency graph (8-connectivity, union–find
  labelling) and its statistics (`bag_stats`, `vertex_scaling`).
* Independent checks, run on the finest basin of each system: box-counting dimension of
  the boundary from the uncertainty exponent (`uncertainty_exponent`) and the merging
  method Wada test (`test_wada_merge`). Per-system results are saved in
  `data/paper_results`.

A full run from scratch takes of the order of half an hour on 8 threads; the magnetic
pendulum dominates.
