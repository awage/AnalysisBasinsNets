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
using Random
using Statistics

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
# 4D system: state = [x, y, vx, vy], N magnets placed on the unit circle
# Parameters: γs — per-magnet strengths, d — height above magnets,
#             α — damping, ω — natural (spring) frequency

struct MagneticPendulum
    magnets::Vector{SVector{2,Float64}}
end

struct MagneticPendulumParams
    γs::Vector{Float64}
    d::Float64
    α::Float64
    ω::Float64
end

@inline function (mp::MagneticPendulum)(du, u, p, t)
    γs, d, α, ω = p.γs, p.d, p.α, p.ω
    x, y, vx, vy = u
    Fx = -ω^2 * x
    Fy = -ω^2 * y
    for (i, mag) in enumerate(mp.magnets)
        mx, my = mag
        r3 = ((x - mx)^2 + (y - my)^2 + d^2)^(3/2)
        Fx += γs[i] * (mx - x) / r3
        Fy += γs[i] * (my - y) / r3
    end
    du[1] = vx
    du[2] = vy
    du[3] = -α * vx + Fx
    du[4] = -α * vy + Fy
    return
end

function magnetic_pendulum(u = [sincos(0.12553*2π)..., 0.0, 0.0];
        γ = 1.0, d = 0.3, α = 0.2, ω = 0.5, N = 3,
        γs = fill(γ, N), diffeq = NamedTuple())
    m = MagneticPendulum([SVector(cos(2π*i/N), sin(2π*i/N)) for i in 1:N])
    p = MagneticPendulumParams(γs, d, α, ω)
    return CoupledODEs(m, u, p; diffeq)
end

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

# ── Relaxed Newton family of Saupe (1987) ────────────────────────────────────
# N_h(z) = z - h p(z)/p'(z), p(z) = z^3 + (λ-1)z - λ, on [-2, 2]^2. SAUPE_DF is
# the box-counting dimension of the Julia set from Saupe's Table IV.

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

function get_basins_saupe(h, N)
    data = cached("saupe", Dict(:h => h, :N => N)) do c
        basins = relaxed_newton_basins(c[:h], c[:N])
        @strdict basins
    end
    return data["basins"]
end
