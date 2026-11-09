# Shared setup for the scripts that reproduce the figures and table of the BAG
# paper. Figures and the table are written directly into ../draft.
using DrWatson
@quickactivate "AnalysisBasinsNets"
using CairoMakie
using Printf
using Statistics

include(srcdir("analyze_nets.jl"))
include(srcdir("paper_systems.jl"))

const DRAFTDIR = joinpath(dirname(projectdir()), "draft")
const FIGDIR = joinpath(DRAFTDIR, "figures")
const TABDIR = joinpath(DRAFTDIR, "tables")
mkpath(FIGDIR); mkpath(TABDIR)

# Resolutions of the dynamical systems. NS is used for the expensive flows
# (magnetic pendulum); NS_LARGE extends it where the basins are cheap enough.
# Only N ≥ NFIT_MIN enters the fits: below it |V| is a few tens of vertices and
# the local slopes fluctuate wildly.
const NS = [50, 75, 100, 150, 200, 300, 400, 500]
const NS_LARGE = vcat(NS, [700, 1000])
const NFIT_MIN = 100

# ── Analysis ──────────────────────────────────────────────────────────────────

"""
    analyze_system(name, getter, Ns; fit_min = NFIT_MIN, checks = true)

Compute the BAG statistics of `getter(N)` for every `N` in `Ns` and fit the
vertex scaling |V| ~ N^α over the resolutions N ≥ `fit_min`. As a reference,
the box-counting dimension of the boundary is obtained on the same grids from
the number of boundary pixels, N_b ~ N^d (the pixel size is ∝ 1/N), and the
uncertainty exponent is α_u = 2 - d. On the finest basin the merging-method
Wada test is run. Results are saved in data/paper_results/<name>.jld2.
"""
function analyze_system(name, getter, Ns; fit_min = NFIT_MIN, checks = true)
    stats = map(Ns) do N
        B = getter(N)
        s = bag_stats(B)
        nb = boundary_pixels(B)
        @printf("%-22s N=%5d  |V|=%7d  <k>=%.3f  NA=%d  wada=%.3f  Nb=%8d  unclassified=%.2e\n",
                name, N, s.nv, s.k, s.NA, s.wada_frac, nb, s.frac_unclassified)
        merge(s, (nb = nb,))
    end
    fit = findall(≥(fit_min), Ns)
    V = [s.nv for s in stats]
    α, b, R2 = vertex_scaling(Ns[fit], V[fit])
    @printf("%-22s α = %.3f  (R² = %.4f, N ≥ %d)\n", name, α, R2, fit_min)
    res = Dict{String,Any}("name" => name, "Ns" => collect(Ns), "stats" => stats,
                           "fit" => fit, "alpha" => α, "b" => b, "R2" => R2)
    Nb = [s.nb for s in stats]
    if all(>(0), Nb[fit])
        d, _, R2d = vertex_scaling(Ns[fit], Nb[fit])
        res["d_box"], res["alpha_u"], res["R2_box"] = d, 2 - d, R2d
    else   # a single basin has no boundary
        res["d_box"], res["alpha_u"], res["R2_box"] = NaN, NaN, NaN
    end
    @printf("%-22s d (boundary box counting) = %.3f  (α_u = %.3f)\n", name, res["d_box"], res["alpha_u"])
    if checks && stats[end].NA ≥ 3
        res["wada_merge"] = wada_merge(getter(last(Ns)))
        @printf("%-22s merging-method Wada test = %s\n", name, res["wada_merge"])
    end
    wsave(datadir("paper_results", "$(name).jld2"), res)
    return res
end

# Number of pixels with an 8-neighbor in a different basin. Counting them at
# resolution N is box counting of the boundary with boxes of one pixel.
function boundary_pixels(B)
    B = fill_unclassified(B)
    return count(CartesianIndices(B)) do I
        any(B[J] != B[I] for J in neighbor_indices(I, size(B)))
    end
end

# Merging method (Daza et al. 2018): fraction of boundary points that fail to
# be on the boundary of all basins after fattening by r pixels. Wada ⇔ → 0.
# test_wada_merge needs labels 1..NA, so labels are compacted first.
function wada_merge(B; r = 1)
    B = fill_unclassified(B)
    ids = sort(unique(B))
    lut = Dict(id => k for (k, id) in enumerate(ids))
    return test_wada_merge(map(b -> lut[b], B), r)
end

# The independent checks only concern the basins of the attractors, so
# unclassified pixels (label ≤ 0) take the label of a classified neighbor.
function fill_unclassified(B)
    B = copy(B)
    while any(<=(0), B)
        todo = findall(<=(0), B)
        for I in todo
            for J in neighbor_indices(I, size(B))
                if B[J] > 0
                    B[I] = B[J]
                    break
                end
            end
        end
        all(B[I] <= 0 for I in todo) && break   # nothing classified at all
    end
    return B
end

load_results(name) = wload(datadir("paper_results", "$(name).jld2"))

# ── Plotting ──────────────────────────────────────────────────────────────────

const BASIN_COLORS = Dict(
    1 => colorant"#d62728", 2 => colorant"#2ca02c", 3 => colorant"#1f4fd6",
    4 => colorant"#ff7f0e")
const UNCLASSIFIED_COLOR = colorant"gray70"

basin_image(B; colors = BASIN_COLORS) =
    [get(colors, b, UNCLASSIFIED_COLOR) for b in B]

function set_theme_paper!()
    set_theme!(Theme(fontsize = 11, figure_padding = 4,
                     Axis = (xgridvisible = false, ygridvisible = false)))
end

function plot_basins(file, xs, ys, B; xlabel = L"x", ylabel = L"y", title = "",
                     colors = BASIN_COLORS, size = (300, 280), extra! = ax -> nothing)
    fig = Figure(; size)
    ax = Axis(fig[1, 1]; xlabel, ylabel, title, aspect = DataAspect())
    image!(ax, first(xs)..last(xs), first(ys)..last(ys), basin_image(B; colors);
           interpolate = false)
    extra!(ax)
    limits!(ax, first(xs), last(xs), first(ys), last(ys))
    save(joinpath(FIGDIR, file), fig; pt_per_unit = 1)
    return fig
end

function plot_scaling(file, series; title = "", size = (320, 280), legend = :lt)
    fig = Figure(; size)
    ax = Axis(fig[1, 1]; xlabel = L"\log N", ylabel = L"\log |V|", title)
    for (k, r) in enumerate(series)
        x = log.(r["Ns"]); y = log.([s.nv for s in r["stats"]])
        fit = get(r, "fit", eachindex(x))
        out = setdiff(eachindex(x), fit)
        c = Makie.wong_colors()[k]
        label = haskey(r, "label") ? @sprintf("%s, α = %.3f", r["label"], r["alpha"]) :
                                     @sprintf("α = %.3f", r["alpha"])
        scatter!(ax, x[fit], y[fit]; color = c, label)
        isempty(out) || scatter!(ax, x[out], y[out]; color = :white, strokecolor = c,
                                 strokewidth = 1)
        lines!(ax, x[fit], r["b"] .+ r["alpha"] .* x[fit]; color = c, label)
    end
    axislegend(ax; position = legend, framevisible = false, merge = true)
    save(joinpath(FIGDIR, file), fig; pt_per_unit = 1)
    return fig
end

set_theme_paper!()
