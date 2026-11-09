# Patch-level description of a categorical map (basins of attraction or any
# land-cover-like raster) through its basin adjacency graph (BAG), and the
# metrics discussed in draft/patchiness_note.tex:
#   first order   α = 2K (Korčák) from the vertex scaling, or 2K on a single grid
#   second order  β (k ~ A^β) = p·γ, degree tail τ ≈ α/(2β)
#   third order   triangles (triple junctions), Wada fraction
#   weighted      bottleneck area a* of the cheapest crossing of the domain
# Requires analyze_nets.jl (label_components, scan_basins, neighbor_indices, bag_stats).

using Graphs
using Random
using Statistics

ols_slope(x, y) = (lx = log.(x); ly = log.(y); cov(lx, ly) / var(lx))

"""
    PatchMap(B)

Patches (8-connected components of equal label) of the map `B`, with their
basin label, area and number of boundary pixels (pixels with an 8-neighbour of
another label) in pixels, centroid, and the BAG `g`.
"""
struct PatchMap
    label::Vector{Int}
    area::Vector{Int}
    perim::Vector{Int}
    cx::Vector{Float64}
    cy::Vector{Float64}
    g::SimpleGraph{Int}
    npix::Int
    dims::Tuple{Int,Int}
    patch::Matrix{Int}
end

function PatchMap(B::AbstractMatrix; background = nothing)
    v, g, labf = scan_basins(B)
    lab = Int.(labf); labf = nothing; M = nv(g)
    area = zeros(Int, M); perim = zeros(Int, M); cx = zeros(M); cy = zeros(M)
    for I in CartesianIndices(B)
        l = lab[I]
        area[l] += 1; cx[l] += I[1]; cy[l] += I[2]
        any(B[J] != B[I] for J in neighbor_indices(I, size(B))) && (perim[l] += 1)
    end
    cx ./= area; cy ./= area
    v = Int.(v)
    background === nothing && return PatchMap(v, area, perim, cx, cy, g, length(B), size(B), lab)
    # Drop the patches of the background label (e.g. the non-forest matrix of a
    # habitat map): the graph keeps the patches of the other classes and their
    # mutual contacts; background pixels get patch index 0.
    keep = findall(!=(background), v)
    g2, vmap = induced_subgraph(g, keep)
    newidx = zeros(Int, M); newidx[vmap] .= 1:length(vmap)
    for i in eachindex(lab); lab[i] = newidx[lab[i]]; end
    PatchMap(v[vmap], area[vmap], perim[vmap], cx[vmap], cy[vmap], g2, length(B), size(B), lab)
end

boundary_pixel_count(pm::PatchMap) = sum(pm.perim)

# ── First order: Korčák law on a single grid ─────────────────────────────────

"""
    korcak_exponent(pm; amin = 10, amax_frac = 1e-3) -> (K, nthresholds)

Korčák exponent of N(A > a) ~ a^(-K) with the thresholds of Imre et al. (2012),
a_i = A_i - 1 over the distinct patch areas A_i in [amin, amax_frac·npix] pixels.
D_K = 2K is the single-grid counterpart of the vertex-scaling exponent α.
"""
function korcak_exponent(pm::PatchMap; amin = 10, amax_frac = 1e-3)
    A = pm.area
    us = sort(unique(A)); us = us[(us .≥ amin) .& (us .≤ amax_frac * pm.npix)]
    length(us) < 3 && return NaN, length(us)
    x = log.(us .- 0.5); y = log.([count(>(u - 1), A) for u in us])
    return -cov(x, y) / var(x), length(us)
end

# ── Second order: degree, area and perimeter ─────────────────────────────────

"""
    contact_exponents(pm; amin = 4) -> (β, p, γ, n)

Log-log least-squares exponents over the patches with area ≥ `amin` pixels and
at least one neighbour: k ~ A^β (degree–area), P ~ A^p (boundary pixels–area,
p = 1/2 for compact patches, 1 for patches made only of boundary) and
k ~ P^γ (neighbours per boundary length). β ≈ p·γ.
"""
function contact_exponents(pm::PatchMap; amin = 4)
    k = degree(pm.g)
    s = (pm.area .≥ amin) .& (pm.perim .> 0) .& (k .> 0)
    n = count(s)
    n < 10 && return NaN, NaN, NaN, n
    A = Float64.(pm.area[s]); P = Float64.(pm.perim[s]); K = Float64.(k[s])
    return ols_slope(A, K), ols_slope(A, P), ols_slope(P, K), n
end

"""
    degree_tail(pm; kmin = 6)

Maximum-likelihood exponent τ of the degree CCDF, P(k' ≥ k) ~ k^(-τ) for k ≥ kmin
(continuous approximation with the usual -1/2 correction).
"""
function degree_tail(pm::PatchMap; kmin = 6)
    k = degree(pm.g); t = k[k .≥ kmin]
    length(t) < 20 && return NaN
    return length(t) / sum(log.(t ./ (kmin - 0.5)))
end

# ── Third order: triple junctions ────────────────────────────────────────────

# Edges join patches of different basins, so every triangle of the BAG is a
# triple contact of three patches of three different basins.
triangle_count(pm::PatchMap) = sum(triangles(pm.g)) ÷ 3

# ── Weighted graph: bottleneck crossing ──────────────────────────────────────

"""
    bottleneck(pm) -> (a_lr, a_tb)

Patches are added by increasing area (union–find on the BAG) until the added
set connects the left and right edges of the map (resp. top and bottom).
The area of the patch that closes the crossing, as a fraction of the map, is
the bottleneck a*: the largest patch that every crossing must traverse.
a* → 0 as the resolution grows when the small patches form a spanning web;
a* stays finite when large patches are barriers.
"""
function bottleneck(pm::PatchMap)
    M = nv(pm.g); N1, N2 = pm.dims; lab = pm.patch
    L = falses(M); R = falses(M); T = falses(M); D = falses(M)
    for j in 1:N2
        lab[1, j] > 0 && (L[lab[1, j]] = true); lab[N1, j] > 0 && (R[lab[N1, j]] = true)
    end
    for i in 1:N1
        lab[i, 1] > 0 && (T[lab[i, 1]] = true); lab[i, N2] > 0 && (D[lab[i, N2]] = true)
    end
    order = sortperm(pm.area)
    out = map(((L, R), (T, D))) do (P, Q)
        parent = collect(1:M); hasP = copy(P); hasQ = copy(Q); added = falses(M)
        for i in order
            added[i] = true
            for j in neighbors(pm.g, i)
                added[j] || continue
                ri, rj = uf_find!(parent, i), uf_find!(parent, j)
                ri == rj && continue
                parent[rj] = ri; hasP[ri] |= hasP[rj]; hasQ[ri] |= hasQ[rj]
            end
            r = uf_find!(parent, i)
            hasP[r] && hasQ[r] && return pm.area[i] / pm.npix
        end
        NaN
    end
    return out[1], out[2]
end

# Wada fraction: share of boundary patches (with at least one neighbour) whose
# neighbours carry every other class. Labels ≤ 0 (unclassified) are not classes.
function wada_fraction(pm::PatchMap)
    NA = count(>(0), unique(pm.label))
    bnd = [i for i in 1:nv(pm.g) if degree(pm.g, i) > 0 && pm.label[i] > 0]
    isempty(bnd) && return NaN, NA
    w = [length(unique(filter(>(0), pm.label[neighbors(pm.g, i)]))) for i in bnd]
    return count(==(NA - 1), w) / length(bnd), NA
end

# Mean Wada neighbour count w̄: average number of other classes adjacent to a
# boundary patch (w̄ = N_A - 1 when every boundary patch touches all classes).
function wada_mean(pm::PatchMap)
    bnd = [i for i in 1:nv(pm.g) if degree(pm.g, i) > 0 && pm.label[i] > 0]
    isempty(bnd) && return NaN
    mean(length(unique(filter(>(0), pm.label[neighbors(pm.g, i)]))) for i in bnd)
end

# ── Distances ────────────────────────────────────────────────────────────────

# Mean hop distance from `nsrc` random patches with at least one neighbour.
function mean_hops(pm::PatchMap; nsrc = 50, rng = Xoshiro(1))
    k = degree(pm.g); cand = findall(>(0), k)
    isempty(cand) && return NaN
    srcs = rand(rng, cand, min(nsrc, length(cand)))
    mean(mean(filter(d -> 0 < d < typemax(Int), gdistances(pm.g, s))) for s in srcs)
end

# ── Collections ──────────────────────────────────────────────────────────────

"""
    fixed_scale_metrics(B) -> NamedTuple

All the metrics that need a single map: number of patches, 2K, β, p, γ, τ,
α/(2β) is left to the caller (needs α), triangles per patch, Wada fraction,
largest degree, mean hop distance, bottleneck a* in both directions.
"""
function fixed_scale_metrics(B; background = nothing)
    pm = PatchMap(B; background)
    K, nK = korcak_exponent(pm)
    β, p, γ, nβ = contact_exponents(pm)
    wada, NA = wada_fraction(pm)
    a_lr, a_tb = bottleneck(pm)
    k = degree(pm.g)
    (; nv = nv(pm.g), nb = boundary_pixel_count(pm), twoK = 2K, nK, β, p, γ, nβ,
       τ = degree_tail(pm), tri_per_v = triangle_count(pm) / nv(pm.g),
       wada, NA, kmax = maximum(k), mean_k = mean(k), isolated = mean(k .== 0),
       hops = mean_hops(pm),
       a_lr, a_tb)
end

"""
    resolution_metrics(get, Ns) -> NamedTuple

Fixed-scale metrics of `get(N)` for every N, and the exponents that need the
resolution sweep: α (|V| ~ N^α), d_b (boundary pixels ~ N^d_b), largest degree
and mean hop distance against |V|. Fits use all the given N.
"""
function resolution_metrics(get, Ns; background = nothing)
    ms = [fixed_scale_metrics(get(N); background) for N in Ns]
    nv = [m.nv for m in ms]
    (; Ns = collect(Ns), per_N = ms, α = ols_slope(Ns, nv),
       d_b = ols_slope(Ns, [m.nb for m in ms]),
       kmax_exp = ols_slope(nv, [m.kmax for m in ms]),
       hops_exp = ols_slope(nv, [m.hops for m in ms]),
       last = last(ms))
end

# ── Lacunarity ───────────────────────────────────────────────────────────────

"""
    gliding_box_lacunarity(mask, r)

Λ(r) = ⟨M²⟩/⟨M⟩² of the mass of `mask` in all r×r boxes (Allain & Cloitre 1991).
"""
function gliding_box_lacunarity(mask::AbstractMatrix{Bool}, r)
    F = Float64.(mask); N1, N2 = size(F)
    S = zeros(N1 + 1, N2 + 1); S[2:end, 2:end] = cumsum(cumsum(F, dims = 1), dims = 2)
    M = [S[i+r, j+r] - S[i, j+r] - S[i+r, j] + S[i, j] for i in 1:N1-r+1, j in 1:N2-r+1]
    mean(M .^ 2) / mean(M)^2
end
