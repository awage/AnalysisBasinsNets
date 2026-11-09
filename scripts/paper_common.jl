# Shared setup for the scripts that reproduce the figures and the table of the
# BAG paper. Figures and the table are written directly into ../draft.
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

fmt(x; d = 3) = (x isa Real && isfinite(x)) ? @sprintf("%.*f", d, x) : "---"

# ── Analysis ──────────────────────────────────────────────────────────────────

"""
    analyze_system(name, getter, Ns; merge_test = false)

BAG statistics of `getter(N)` for every `N` in `Ns`, together with the number of
boundary pixels `nb` (box counting of the boundary with one-pixel boxes) and,
on the finest grid only, the single-grid Korčák dimension 2K. With
`merge_test = true` the merging-method Wada test is run on the finest grid.
Exponents are fitted afterwards with `exponents`. Results are saved in
data/paper_results/<name>.jld2.
"""
function analyze_system(name, getter, Ns; merge_test = false)
    stats = map(Ns) do N
        B = getter(N)
        s = merge(bag_stats(B; korcak = N == last(Ns)), (nb = boundary_pixels(B),))
        @printf("%-14s N=%5d  |V|=%7d  <k>=%.2f  kmax=%6d  wada=%.3f  w̄=%.3f  Nb/N²=%.3f  unclassified=%.1e\n",
                name, N, s.nv, s.k, s.kmax, s.wada_frac, s.w, s.nb / N^2, s.frac_unclassified)
        s
    end
    res = Dict{String,Any}("name" => name, "Ns" => collect(Ns), "stats" => stats)
    if merge_test
        res["wada_merge"] = wada_merge(getter(last(Ns)))
        @printf("%-14s merging-method Wada test = %.4f\n", name, res["wada_merge"])
    end
    wsave(datadir("paper_results", "$(name).jld2"), res)
    return res
end

"""
    exponents(res; Nmin = 0) -> NamedTuple

Exponents fitted over the grids N ≥ `Nmin` of an `analyze_system` result:
vertex scaling |V| ~ N^α (with intercept `b` and `R2`), boundary pixels
N_b ~ N^d_b, largest degree kmax ~ |V|^x.
"""
function exponents(res; Nmin = 0)
    Ns = res["Ns"]; fit = findall(≥(Nmin), Ns); s = res["stats"][fit]
    V = [x.nv for x in s]
    α, b, R2 = vertex_scaling(Ns[fit], V)
    d_b, = vertex_scaling(Ns[fit], [x.nb for x in s])
    x_kmax, = vertex_scaling(V, [x.kmax for x in s])
    (; α, b, R2, d_b, x_kmax, fit, Ns = Ns[fit], stats = s)
end

# Number of pixels with an 8-neighbour in a different basin. Counting them at
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

# The boundary-pixel count and the merging test only concern the basins of the
# attractors, so unclassified pixels (label ≤ 0) take the label of a
# classified neighbour.
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

# ── Plotting ──────────────────────────────────────────────────────────────────

const BASIN_COLORS = Dict(
    1 => colorant"#d62728", 2 => colorant"#2ca02c", 3 => colorant"#1f4fd6",
    4 => colorant"#ff7f0e")
const UNCLASSIFIED_COLOR = colorant"gray70"

basin_image(B; colors = BASIN_COLORS) =
    [get(colors, b, UNCLASSIFIED_COLOR) for b in B]

set_theme!(Theme(fontsize = 11, figure_padding = 4,
                 Axis = (xgridvisible = false, ygridvisible = false)))
