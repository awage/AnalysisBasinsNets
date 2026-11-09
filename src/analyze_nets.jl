using Graphs
using Statistics

# Return CartesianIndices of all neighbors of `box` within `dims` bounds
function neighbor_indices(box, dims)
    CartesianIndices(ntuple(i -> max(1, box[i]-1):min(dims[i], box[i]+1), length(dims)))
end

function cart_to_lin(k, dims)
    return k[1] + (k[2]-1)*dims[1]
end

# ── Union-find used to label connected components of the pixel graph ─────────

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

# Label the 8-connected components of equal basin value. Components are
# numbered 1, 2, ... in order of first appearance (column-major scan).
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

# Build the basin adjacency graph (BAG). Returns the basin value of each
# component, the BAG and the component-numbered basin.
function scan_basins(basins)
    dims = size(basins)
    lab, M = label_components(basins)

    basins_values = zeros(M)
    for i in eachindex(lab)
        basins_values[lab[i]] = basins[i]
    end

    # Connect adjacent components (they belong to different basins by construction)
    basins_graph = SimpleGraph(M)
    for box in CartesianIndices(basins)
        a = lab[box]
        for k in neighbor_indices(box, dims)
            b = lab[k]
            a != b && add_edge!(basins_graph, a, b)
        end
    end

    return basins_values, basins_graph, Float64.(lab)
end

function wada_neighbors(g, v)
    wd = zeros(length(v))
    for (k, vtx) in enumerate(g.fadjlist)
        wd[k] = length(unique(v[vtx]))
    end
    return wd
end

function scan_boundary!(G, basins)
    dims = size(basins)
    box_indices = CartesianIndices(basins)
    for box in box_indices
        nb = cart_to_lin(box, dims)
        for k in neighbor_indices(box, dims)
            n = cart_to_lin(k, dims)
            if basins[box] != basins[k] && n != nb
                add_edge!(G, n, nb)
            end
        end
    end
end

# ── Graph statistics used in the paper ────────────────────────────────────────

"""
    bag_stats(basins) -> NamedTuple

Summary statistics of the basin adjacency graph of `basins`:
number of vertices `nv` and edges `ne`, mean degree `k`, number of
basins `NA`, and the Wada fraction: the fraction of boundary components
(components with at least one neighbor) whose neighbors carry all the
`NA - 1` other basin labels. `w` is the mean Wada neighbor count over the
boundary components. Pixels with a label ≤ 0 (initial conditions that did
not converge) are kept as regions of the graph but are not counted as a
basin; `frac_unclassified` is their fraction.
"""
function bag_stats(basins)
    v, g, _ = scan_basins(basins)
    NA = count(>(0), unique(basins))
    deg = degree(g)
    wd = [length(unique(filter(>(0), v[nbrs]))) for nbrs in g.fadjlist]
    bnd = (deg .> 0) .& (v .> 0)
    nb = count(bnd)
    return (
        nv = nv(g),
        ne = ne(g),
        k = mean(deg),
        NA = NA,
        wada_frac = nb > 0 ? count(wd[bnd] .== NA - 1) / nb : NaN,
        w = nb > 0 ? mean(wd[bnd]) : NaN,
        frac_unclassified = count(<=(0), basins) / length(basins),
    )
end

"""
    vertex_scaling(Ns, V) -> (α, b, R²)

Least-squares fit of log|V| = b + α log N.
"""
function vertex_scaling(Ns, V)
    x = log.(Ns); y = log.(V)
    α = cov(x, y) / var(x)
    b = mean(y) - α * mean(x)
    R2 = cor(x, y)^2
    return α, b, R2
end
