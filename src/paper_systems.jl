# Basins used in the BAG paper (draft in ../draft). Every generator returns an
# N×N integer matrix indexed [i_x, j_y] with basin labels 1, 2, ...; negative
# labels mark initial conditions that were not classified.
#
# Basins of the dynamical systems are computed in two steps so that the labels
# are consistent across resolutions: the attractors are located once
# (`*_attractors`), then each N×N grid is classified with a threaded
# AttractorsViaProximity mapper (`proximity_basins`). Results are cached in
# data/paper_basins with DrWatson.produce_or_load.

using DrWatson
using OrdinaryDiffEqVerner
using Attractors
using LinearAlgebra
using Random
using Statistics

include(srcdir("basins_fun.jl"))   # magnetic_pendulum

# ── Sierpiński benchmarks ─────────────────────────────────────────────────────

# Level-m triangle on an N = 2^m grid (Pascal's triangle mod 2). Filled cells
# are colored 1–3 by the last binary digit of their address, so that touching
# sub-triangles have different colors; the holes and background are label 4.
function sierpinski_triangle(m)
    N = 2^m
    B = fill(4, N, N)
    for x in 0:N-1, y in 0:N-1
        if (x & y) == x
            B[x+1, y+1] = 1 + (x & 1) + (y & 1)
        end
    end
    return B
end

# Level-m carpet on an N = 3^m grid: filled = 1, holes = 2.
function sierpinski_carpet(m)
    N = 3^m
    B = ones(Int, N, N)
    for i in 0:N-1, j in 0:N-1
        a, b = i, j
        while a > 0 || b > 0
            if a % 3 == 1 && b % 3 == 1
                B[i+1, j+1] = 2
                break
            end
            a ÷= 3; b ÷= 3
        end
    end
    return B
end

# ── Generic machinery ─────────────────────────────────────────────────────────

# Classify every point of the grid (xs, ys) with a proximity mapper built by
# `make_mapper()` (one per task, since mappers are stateful).
function proximity_basins(make_mapper, xs, ys; to_state = (x, y) -> SVector(x, y))
    B = zeros(Int, length(xs), length(ys))
    chunk = max(1, length(ys) ÷ (4 * Threads.nthreads()))
    @sync for cols in Iterators.partition(eachindex(ys), chunk)
        Threads.@spawn begin
            mapper = make_mapper()
            for j in cols, i in eachindex(xs)
                B[i, j] = mapper(to_state(xs[i], ys[j]))
            end
        end
    end
    return B
end

# Locate attractors with the recurrences algorithm from a coarse set of
# initial conditions.
function recurrence_attractors(ds, rec_grid, ics; kwargs...)
    mapper = AttractorsViaRecurrences(ds, rec_grid; kwargs...)
    for u in ics
        mapper(u)
    end
    return extract_attractors(mapper)
end

function cached(f, name, config)
    data, _ = produce_or_load(config, datadir("paper_basins"); prefix = name) do c
        f(c)
    end
    return data
end

# ── Newton map for z^3 - 1 on [-1, 1]^2 (window of Epureanu & Greenside) ──────────────────────────────────────

function newton_rule(u, p, n)
    z = complex(u[1], u[2])
    abs2(z) < 1e-24 && return SVector(1e15, 1e15)  # z = 0 is sent to infinity
    z = z - (z^3 - 1) / (3z^2)
    return SVector(real(z), imag(z))
end

const NEWTON_BOX = (-1.0, 1.0)

function newton_attractors()
    roots = [SVector(cospi(2k/3), sinpi(2k/3)) for k in 0:2]
    return Dict(k => StateSpaceSet([roots[k]]) for k in 1:3)
end

function get_basins_newton(N)
    data = cached("newton", Dict(:N => N, :L => last(NEWTON_BOX))) do c
        ds = DeterministicIteratedMap(newton_rule, SVector(0.1, 0.2))
        att = newton_attractors()
        xs = ys = range(NEWTON_BOX..., length = c[:N])
        basins = proximity_basins(xs, ys) do
            AttractorsViaProximity(deepcopy(ds), att, 1e-3;
                                   consecutive_lost_steps = 10_000, horizon_limit = 1e14)
        end
        @strdict basins
    end
    return data["basins"]
end

# Root reached by Newton's method from z (1–3), 0 if it does not converge. Used
# for the grid-free uncertainty exponent below.
function newton_root(z::Complex)
    roots = (complex(1.0, 0.0), cis(2π/3), cis(4π/3))
    for _ in 1:500
        abs2(z) < 1e-24 && return 0
        z -= (z^3 - 1) / (3z^2)
        for k in 1:3
            abs2(z - roots[k]) < 1e-6 && return k
        end
    end
    return 0
end

"""
    grebogi_exponent(label, box; εs, M, seed) -> (α_u, εs, f)

Uncertainty exponent computed as in McDonald et al. (1985), without a grid:
`M` random initial conditions z in `box`², a point is uncertain if `label`
differs at z ± ε along the first axis. f(ε) ~ ε^α_u, and d = 2 - α_u.
"""
function grebogi_exponent(label, box; εs = 10.0 .^ (-2:-0.5:-6.5), M = 400_000, seed = 1)
    a, b = box
    f = map(εs) do ε
        c = Threads.Atomic{Int}(0)
        Threads.@threads for i in 1:M
            rng = Xoshiro(hash((seed, i)))
            x0, y0 = a + (b - a) * rand(rng), a + (b - a) * rand(rng)
            l = label(x0, y0)
            (l != label(x0 - ε, y0) || l != label(x0 + ε, y0)) && Threads.atomic_add!(c, 1)
        end
        c[] / M
    end
    lε, lf = log.(εs), log.(f)
    return cov(lε, lf) / var(lε), εs, f
end

# ── Duffing oscillator  x'' + δx' - x + x^3 = F sin(ωt) ──────────────────────

function duffing_sin(u, p, t)
    δ, F, ω = p
    return SVector(u[2], -δ*u[2] + u[1] - u[1]^3 + F*sin(ω*t))
end

# Wada parameters of Aguirre & Sanjuán, Physica D 171, 41 (2002), as used in
# Daza et al., PLoS ONE 13, e0194926 (2018): δ = 0.15, F = 0.245, ω = 1 on
# [-2, 2]^2. Two period-1 attractors and one period-3 attractor, Wada boundary.
# With the sine forcing the stroboscopic section at t = 0 mod 2π reproduces the
# published attractors, P1R ≈ (0.815, 0.242), P1L ≈ (-0.933, 0.299) and
# P3 ≈ (-1.412, -0.137), (-0.354, -0.614), (0.645, -0.464). (PLoS ONE prints
# cos ωt, but its attractors are those of the sine forcing.)
const DUFFING_δ = 0.15
const DUFFING_F = 0.245
const DUFFING_ω = 1.0
# The proximity mapper labels an orbit the first time it comes ε-close to an
# attractor. The default ε (1/10 of the distance between attractors) lets the
# long chaotic transients be captured early by the period-3 orbit; ε = 1e-3
# gives the same labels as 1e-4 on a 120×120 grid.
const DUFFING_PROX_ε = 1e-3
const DUFFING_BOX = (-2.0, 2.0)

duffing_map(F) = StroboscopicMap(
    CoupledODEs(duffing_sin, SVector(0.1, 0.1), [DUFFING_δ, F, DUFFING_ω];
                diffeq = (alg = Vern9(), reltol = 1e-8, abstol = 1e-8)),
    2π / DUFFING_ω)

function duffing_attractors(F)
    data = cached("duffing_attractors", Dict(:F => F, :δ => DUFFING_δ, :ω => DUFFING_ω, :forcing => "sin")) do c
        smap = duffing_map(F)
        rec = (range(-3, 3, length = 601), range(-3, 3, length = 601))
        g = range(DUFFING_BOX..., length = 40)
        ics = [SVector(x, y) for x in g for y in g]
        attractors = recurrence_attractors(smap, rec, ics; consecutive_recurrences = 200)
        @strdict attractors
    end
    return data["attractors"]
end

function get_basins_duffing(F, N)
    data = cached("duffing", Dict(:F => F, :N => N, :δ => DUFFING_δ, :ω => DUFFING_ω, :forcing => "sin",
                           :ε => DUFFING_PROX_ε)) do c
        att = duffing_attractors(F)
        smap = duffing_map(F)
        xs = ys = range(DUFFING_BOX..., length = c[:N])
        basins = proximity_basins(xs, ys) do
            AttractorsViaProximity(deepcopy(smap), att; ε = DUFFING_PROX_ε,
                                   consecutive_lost_steps = 20_000)
        end
        @strdict basins
    end
    return data["basins"]
end

# ── Magnetic pendulum, three magnets on the unit circle ──────────────────────

const MAG_PARAMS = (γ = 1.0, d = 0.3, α = 0.05, ω = 0.5, N = 3)
const MAG_BOX = (-2.0, 2.0)

mag_system() = magnetic_pendulum(; MAG_PARAMS...,
    diffeq = (alg = Vern9(), reltol = 1e-8, abstol = 1e-8))

mag_magnets() = [SVector(cospi(2i/3), sinpi(2i/3)) for i in 1:MAG_PARAMS.N]

# The attractors are the stable equilibria next to each magnet: integrate from
# rest on top of each magnet until the pendulum settles.
function mag_attractors()
    ds = mag_system()
    att = Dict{Int,StateSpaceSet}()
    for (k, m) in enumerate(mag_magnets())
        reinit!(ds, [m..., 0.0, 0.0])
        step!(ds, 2000.0)
        att[k] = StateSpaceSet([SVector{4}(current_state(ds))])
    end
    return att
end

function get_basins_magnetic(N)
    data = cached("magnetic", merge(Dict{Symbol,Any}(pairs(MAG_PARAMS)), Dict(:res => N))) do c
        att = mag_attractors()
        ds = mag_system()
        xs = ys = range(MAG_BOX..., length = c[:res])
        basins = proximity_basins(xs, ys; to_state = (x, y) -> [x, y, 0.0, 0.0]) do
            AttractorsViaProximity(deepcopy(ds), att; Δt = 0.1,
                                   consecutive_lost_steps = 50_000)
        end
        @strdict basins
    end
    return data["basins"]
end

# ── Quasiperiodically forced map (Feudel, Witt, Lai, Grebogi 1998) ───────────
# State (x, u); M = f∘f∘f with f(x) = r x (1 - x).

function feudel_rule(z, p, n)
    a, r = p
    f(x) = r * x * (1 - x)
    x, u = z
    ν = (sqrt(5.0) - 1) / 2
    return SVector(f(f(f(x))) + a * cospi(2u), mod(u + ν, 1.0))
end

const FEUDEL_PARAMS = (a = 0.0015, r = 3.833)
# As for the Duffing oscillator, the default proximity ε captures transients
# early (64% agreement with ε = 1e-3 on a 150×150 grid); 1e-3 and 1e-4 agree.
const FEUDEL_PROX_ε = 1e-3

feudel_map() = DeterministicIteratedMap(feudel_rule, SVector(0.5, 0.0),
                                        [FEUDEL_PARAMS.a, FEUDEL_PARAMS.r])

function feudel_attractors()
    data = cached("feudel_attractors", Dict(pairs(FEUDEL_PARAMS))) do c
        rec = (range(-0.1, 1.1, length = 1201), range(0, 1, length = 1001))
        g = range(0, 1, length = 30)
        ics = [SVector(x, u) for x in g for u in g]
        attractors = recurrence_attractors(feudel_map(), rec, ics;
                                           consecutive_recurrences = 2000)
        @strdict attractors
    end
    return data["attractors"]
end

# Basins on the (u, x) plane: rows are u (phase), columns are x (state).
function get_basins_feudel(N)
    data = cached("feudel", merge(Dict{Symbol,Any}(pairs(FEUDEL_PARAMS)),
                                  Dict(:N => N, :ε => FEUDEL_PROX_ε))) do c
        att = feudel_attractors()
        ds = feudel_map()
        us = xs = range(0, 1, length = c[:N])
        basins = proximity_basins(us, xs; to_state = (u, x) -> SVector(x, u)) do
            AttractorsViaProximity(deepcopy(ds), att; ε = FEUDEL_PROX_ε,
                                   consecutive_lost_steps = 100_000)
        end
        @strdict basins
    end
    return data["basins"]
end
