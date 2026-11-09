# Reference mosaics for draft/patchiness_note.tex: a neutral random mosaic with
# tunable fragmentation, two carpets with equal dimension and different
# lacunarity, the relaxed Newton family of Saupe (1987) and quadratic Julia sets.

using Random

# ── Neutral random quadtree mosaic ───────────────────────────────────────────
# A cell at depth j splits into four with probability q (the first two levels
# always split). Split decisions and colours are hashed from the cell address,
# so the same mosaic is rasterised at every resolution 2^d. With colors = 0
# every leaf gets its own label (no merging); otherwise leaves are coloured at
# random with `colors` labels and equal neighbours merge into irregular patches.
# Expected scaling: α = log2(4q); for unique labels β = max(0, log2(2q)/2).

quadtree_splits(q, j, a, b, seed) = j < 2 || rand(Xoshiro(hash((seed, j, a, b)))) < q

function quadtree_mosaic(q, d; colors = 3, seed = 7)
    N = 2^d; B = zeros(Int, N, N); id = Ref(0)
    function fill_cell!(j, a, b)
        if j < d && quadtree_splits(q, j, a, b, seed)
            for (da, db) in ((0, 0), (1, 0), (0, 1), (1, 1))
                fill_cell!(j + 1, 2a + da, 2b + db)
            end
        else
            s = 2^(d - j); id[] += 1
            c = colors == 0 ? id[] : 1 + Int(hash((seed, j, a, b, :colour)) % colors)
            B[a*s+1:(a+1)*s, b*s+1:(b+1)*s] .= c
        end
    end
    fill_cell!(0, 0, 0)
    return B
end

# ── Carpets with equal dimension, different lacunarity ───────────────────────
# 7×7 generators keeping 40 of 49 cells (D = log 40 / log 7 ≈ 1.896):
# :block removes the central 3×3 square (one large hole per tile),
# :dust removes 9 isolated interior cells (nine small holes per tile).

const CARPET_GENERATORS = Dict(
    :block => [!(2 ≤ i ≤ 4 && 2 ≤ j ≤ 4) for i in 0:6, j in 0:6],
    :dust  => [!(i in (1, 3, 5) && j in (1, 3, 5)) for i in 0:6, j in 0:6],
)

# Level-m carpet on a 7^m grid: filled = 2, holes = 1.
function carpet7(kind::Symbol, m)
    gen = CARPET_GENERATORS[kind]
    B = trues(1, 1)
    for _ in 1:m
        n = size(B, 1); C = falses(7n, 7n)
        for I in CartesianIndices(gen), J in CartesianIndices(B)
            C[(I[1]-1)*n + J[1], (I[2]-1)*n + J[2]] = gen[I] && B[J]
        end
        B = C
    end
    return Int.(B) .+ 1
end

# ── Relaxed Newton family (Saupe 1987, Table IV) ─────────────────────────────
# N_h(z) = z - h p(z)/p'(z), p(z) = z^3 + (λ-1)z - λ, λ = 1.1 - 0.1i, on [-2, 2]^2.
# Box-counting dimension of the Julia set from Saupe's Table IV:

const SAUPE_λ = 1.1 - 0.1im
const SAUPE_ROOTS = let s = sqrt(1 - 4SAUPE_λ); (complex(1.0), (-1 + s) / 2, (-1 - s) / 2) end
const SAUPE_DF = Dict(0.1 => 1.242, 0.2 => 1.254, 0.3 => 1.290, 0.4 => 1.309, 0.5 => 1.326,
                      0.6 => 1.345, 0.7 => 1.360, 0.8 => 1.369, 0.9 => 1.401, 1.0 => 1.424,
                      1.1 => 1.443, 1.2 => 1.464, 1.3 => 1.485, 1.4 => 1.512, 1.5 => 1.536,
                      1.6 => 1.565, 1.7 => 1.603, 1.8 => 1.638, 1.9 => 1.678)

function relaxed_newton_root(z, h; maxit = 20_000)
    λ = SAUPE_λ
    for _ in 1:maxit
        dp = 3z^2 + (λ - 1)
        abs2(dp) < 1e-30 && return 0
        z -= h * (z^3 + (λ - 1) * z - λ) / dp
        for k in 1:3
            abs2(z - SAUPE_ROOTS[k]) < 1e-12 && return k
        end
    end
    return 0
end

function relaxed_newton_basins(h, N; L = 2.0)
    r = range(-L, L, length = N); B = zeros(Int, N, N)
    Threads.@threads for j in 1:N
        for i in 1:N
            B[i, j] = relaxed_newton_root(complex(r[i], r[j]), h)
        end
    end
    return B
end

# ── Quadratic Julia sets z^2 + c ─────────────────────────────────────────────
# Label 1: escapes, 2: bounded. The basin of ∞ is connected, so these maps test
# the case where the boundary is rough but patches do not proliferate.

function julia_basins(c, N; L = 1.6, maxit = 5000)
    r = range(-L, L, length = N); B = zeros(Int, N, N)
    Threads.@threads for j in 1:N
        for i in 1:N
            z = complex(r[i], r[j]); l = 2
            for _ in 1:maxit
                z = z * z + c
                abs2(z) > 4 && (l = 1; break)
            end
            B[i, j] = l
        end
    end
    return B
end

const JULIA_CASES = [   # (name, c, box-counting dimension from the literature or "")
    ("Douady rabbit", -0.12 + 0.74im, "1.39"),
    ("basilica",      -1.0 + 0.0im,   "1.27"),
    ("quasicircle",   -0.5 + 0.3im,   ""),
]
