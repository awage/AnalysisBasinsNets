# Basin adjacency graph (BAG) of a basin diagram and the statistics used in the
# paper (../draft/main.tex): vertex count, degrees, Wada averages and the
# single-grid Korčák exponent.
using Graphs
using Statistics

# CartesianIndices of the 8-neighbourhood of `box` (itself included) within `dims`.
function neighbor_indices(box, dims)
    CartesianIndices(ntuple(i -> max(1, box[i]-1):min(dims[i], box[i]+1), length(dims)))
end

cart_to_lin(k, dims) = k[1] + (k[2]-1)*dims[1]

# ── Union–find labelling of the 8-connected components ───────────────────────

function uf_find!(parent, i)
    while parent[i] != i
        parent[i] = parent[parent[i]]
        i = parent[i]
    end
    return i
end

function uf_union!(parent, i, j)
    ri, rj = uf_find!(parent, i), uf_find!(parent, j)
    ri == rj && return
    ri < rj ? (parent[rj] = ri) : (parent[ri] = rj)
end

# Components are numbered 1, 2, ... in order of first appearance (column-major).
function label_components(basins)
    dims = size(basins)
    parent = collect(1:length(basins))
    for box in CartesianIndices(basins)
        nb = cart_to_lin(box, dims)
        for k in neighbor_indices(box, dims)
            n = cart_to_lin(k, dims)
            if n < nb && basins[box] == basins[k]
                uf_union!(parent, n, nb)
            end
        end
    end
    lab = zeros(Int, dims)
    root_to_lab = Dict{Int,Int}()
    for i in eachindex(lab)
        r = uf_find!(parent, i)
        lab[i] = get!(root_to_lab, r, length(root_to_lab) + 1)
    end
    return lab, length(root_to_lab)
end

"""
    scan_basins(basins) -> (label, g, lab)

BAG of `basins`: basin label of each component, the graph `g` (edges join
8-adjacent components, which belong to different basins by construction) and
the matrix of component indices.
"""
function scan_basins(basins)
    dims = size(basins)
    lab, M = label_components(basins)
    label = zeros(Int, M)
    for i in eachindex(lab)
        label[lab[i]] = basins[i]
    end
    g = SimpleGraph(M)
    for box in CartesianIndices(basins)
        a = lab[box]
        for k in neighbor_indices(box, dims)
            b = lab[k]
            a != b && add_edge!(g, a, b)
        end
    end
    return label, g, lab
end

# ── Statistics ───────────────────────────────────────────────────────────────

"""
    korcak_exponent(area, npix; amin = 10, amax_frac = 1e-3) -> (K, nthresholds)

Korčák exponent of N(A > a) ~ a^(-K) with the thresholds of Imre et al. (2012),
a_i = A_i - 1 over the distinct component areas A_i in [amin, amax_frac·npix]
pixels. 2K is the single-grid counterpart of the vertex-scaling exponent α.
"""
function korcak_exponent(area, npix; amin = 10, amax_frac = 1e-3)
    us = sort(unique(area)); us = us[(us .≥ amin) .& (us .≤ amax_frac * npix)]
    length(us) < 3 && return NaN, length(us)
    x = log.(us .- 0.5); y = log.([count(>(u - 1), area) for u in us])
    return -cov(x, y) / var(x), length(us)
end

# Number of pixels of each component.
function component_areas(lab, M)
    a = zeros(Int, M)
    for l in lab
        a[l] += 1
    end
    return a
end

"""
    bag_stats(basins; korcak = false) -> NamedTuple

Statistics of the BAG of `basins`: number of vertices `nv` and edges `ne`,
mean and largest degree `k`, `kmax`, number of basins `NA`, Wada fraction
(share of boundary components, i.e. with at least one neighbour, adjacent to
all `NA - 1` other basins) and mean basin-neighbour count `w` over the
boundary components. Pixels with a label ≤ 0 (not classified) are kept as
regions of the graph but are not a basin. With `korcak = true` the
single-grid Korčák dimension `twoK` = 2K is also computed.
"""
function bag_stats(basins; korcak = false)
    v, g, lab = scan_basins(basins)
    NA = count(>(0), unique(basins))
    deg = degree(g)
    wd = [length(unique(filter(>(0), v[nbrs]))) for nbrs in g.fadjlist]
    bnd = (deg .> 0) .& (v .> 0)
    nb = count(bnd)
    twoK = korcak ? 2 * first(korcak_exponent(component_areas(lab, nv(g)), length(basins))) : NaN
    return (
        nv = nv(g),
        ne = ne(g),
        k = mean(deg),
        kmax = maximum(deg),
        NA = NA,
        wada_frac = nb > 0 ? count(wd[bnd] .== NA - 1) / nb : NaN,
        w = nb > 0 ? mean(wd[bnd]) : NaN,
        twoK = twoK,
        frac_unclassified = count(<=(0), basins) / length(basins),
    )
end

"""
    vertex_scaling(Ns, V) -> (α, b, R²)

Least-squares fit of log V = b + α log N.
"""
function vertex_scaling(Ns, V)
    x = log.(Ns); y = log.(V)
    α = cov(x, y) / var(x)
    b = mean(y) - α * mean(x)
    R2 = cor(x, y)^2
    return α, b, R2
end
