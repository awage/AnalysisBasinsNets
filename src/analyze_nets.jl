using Graphs

# Return CartesianIndices of all neighbors of `box` within `dims` bounds
function neighbor_indices(box, dims)
    CartesianIndices(ntuple(i -> max(1, box[i]-1):min(dims[i], box[i]+1), length(dims)))
end

function cart_to_lin(k, dims)
    return k[1] + (k[2]-1)*dims[1]
end

function scan_basins(basins)
    dims = size(basins)
    G = SimpleGraph()
    add_vertices!(G, length(basins))
    box_indices = CartesianIndices(basins)

    # Connect pixels belonging to the same basin
    for box in box_indices
        nb = cart_to_lin(box, dims)
        for k in neighbor_indices(box, dims)
            n = cart_to_lin(k, dims)
            if basins[box] == basins[k] && n != nb
                add_edge!(G, n, nb)
            end
        end
    end

    # Compute connected components and label each region
    components = connected_components(G)
    basins_numbered = zeros(size(basins))
    basins_values = zeros(length(components))
    basins_graph = SimpleGraph(length(components))

    for (k, cmpt) in enumerate(components)
        basins_values[k] = basins[cmpt[1]]
        for vtx in cmpt
            basins_numbered[vtx] = k
        end
    end

    # Connect adjacent components that belong to different basins
    for box in box_indices
        nb = basins_numbered[cart_to_lin(box, dims)]
        for k in neighbor_indices(box, dims)
            n = basins_numbered[cart_to_lin(k, dims)]
            if basins_numbered[box] != basins_numbered[k] && n != nb
                add_edge!(basins_graph, n, nb)
            end
        end
    end

    return basins_values, basins_graph, basins_numbered
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
