using DrWatson
@quickactivate 
using Statistics
using Graphs


include(srcdir("analyze_nets.jl"))
include(srcdir("basins_fun.jl"))


# Compute and print summary statistics for a given basin array
function analyze_basin(basins)
    v, g, b = scan_basins(basins)
    md, sd, nd = mean(degree(g)), std(degree(g)), length(v)
    w  = mean(wada_neighbors(g, v))
    cc = mean(local_clustering_coefficient(g))
    @show md, sd, nd
    @show w
    @show cc
    return v, g, b
end

# See how graph properties grow with boundary resolution → fractal dimension
function compute_stats()
    F = 0.128; ω = 1.106
    res_v = range(100, 500, step=100)
    N = length(res_v)
    md = zeros(N); sd = zeros(N); nd = zeros(N); wd = zeros(N)
    for (k, res) in enumerate(res_v)
        basins = get_basins_duf(F, ω, res)
        v, g, b = scan_basins(basins)
        @show md[k], sd[k], nd[k] = mean(degree(g)), std(degree(g)), length(v)
        @show wd[k] = mean(wada_neighbors(g, v))
    end
    return res_v, md, sd, nd, wd
end


# ── Duffing oscillator cases ──────────────────────────────────────────────────

# Smooth boundary, 2 attractors
println("=== Duffing: smooth boundary, 2 attractors (F=0.1, ω=0.2) ===")
analyze_basin(get_basins_duf(0.1, 0.2, 500))

# Fractal boundary, 2 attractors
println("=== Duffing: fractal boundary, 2 attractors (F=0.2, ω=1.0) ===")
analyze_basin(get_basins_duf(0.2, 1.0, 500))

# Fractal boundary, 4 attractors (Wada candidate)
println("=== Duffing: fractal boundary, 4 attractors (F=0.128, ω=1.106) ===")
analyze_basin(get_basins_duf(0.128, 1.106, 500))

# Fractal boundary, 4 attractors (large forcing)
println("=== Duffing: fractal boundary, 4 attractors (F=2.0, ω=1.0) ===")
analyze_basin(get_basins_duf(2.0, 1.0, 500))

# ── Other systems ─────────────────────────────────────────────────────────────

# Fractal Wada boundary, 3 attractors — quasiperiodically forced map (Feudel 1998)
println("=== Feudel map: fractal Wada boundary, 3 attractors ===")
analyze_basin(get_basins_map(500))

# Fractal boundary, 4 attractors — magnetic pendulum
println("=== Magnetic pendulum: fractal boundary, 4 attractors ===")
analyze_basin(get_bas_mag_pend(500))

# Fractal boundary, 3 attractors — Newton's method roots
println("=== Newton map: fractal boundary, 3 attractors ===")
analyze_basin(get_basins_newt_map(500))
