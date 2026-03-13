using OrdinaryDiffEq
using Attractors


# ── Duffing oscillator ────────────────────────────────────────────────────────

@inline @inbounds function duffing(u, p, t)
    d = p[1]; F = p[2]; omega = p[3]
    du1 = u[2]
    du2 = -d*u[2] + u[1] - u[1]^3 + F*sin(omega*t)
    return SVector{2}(du1, du2)
end

function _get_basins_duf(d)
    @unpack F, ω, res = d
    ds   = CoupledODEs(duffing, rand(2), [0.15, F, ω]; diffeq=(reltol=1e-8, alg=Vern9()))
    smap = StroboscopicMap(ds, 2π/ω)
    xg   = yg = range(-3, 3, length=res)
    mapper = AttractorsViaRecurrences(smap, (xg, yg))
    basins, attractors = basins_of_attraction(mapper, (xg, yg))
    return @strdict basins
end

function get_basins_duf(F, ω, res)
    data, file = produce_or_load(
        datadir("basins"),
        @dict(F, ω, res),
        _get_basins_duf;
        prefix = "basin_duf"
    )
    @unpack basins = data
    return basins
end


# ── Quasiperiodically forced chaotic map ──────────────────────────────────────
# Basin bifurcation in quasiperiodically forced systems
# Feudel, Witt, Lai, Grebogi — PRE 58, 1998

function chaotic_map!(dz, z, p, n)
    xn = z[1]
    θ  = z[2]
    a  = p[1]
    ω  = (sqrt(5.0) - 1.0) / 2.0
    r  = p[2]
    f(x) = r * x * (1.0 - x)
    M = reduce(∘, fill(f, 3))
    dz[1] = M(xn) + a * cos(2π * θ)
    dz[2] = mod(θ + ω, 1.0)
    return
end

function get_basins_map(res)
    ds     = DeterministicIteratedMap(chaotic_map!, [1.0, 0.0], [0.0015, 3.833])
    θg     = range(0.0, 1.0, length=res)
    xg     = range(0.0, 1.0, length=res)
    mapper = AttractorsViaRecurrences(ds, (θg, xg))
    basins, attractors = basins_of_attraction(mapper, (θg, xg))
    return basins
end


# ── Magnetic pendulum ─────────────────────────────────────────────────────────
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

function _get_bas_mag_pend(d)
    @unpack res = d
    diffeq = (alg=Vern9(), reltol=1e-6, maxiters=Int(1e8))
    ds     = magnetic_pendulum(γ=1.0, d=0.3, α=0.2, ω=0.5, N=4; diffeq)
    xg_rec = yg_rec = range(-10.0, 10.0, length=5001)
    psys   = ProjectedDynamicalSystem(ds, [1, 2], [0.0, 0.0])
    mapper = AttractorsViaRecurrences(psys, (xg_rec, yg_rec); Δt=0.1)
    xg     = yg = range(-2.0, 2.0, length=res)
    basins, attractors = basins_of_attraction(mapper, (xg, yg))
    return @strdict basins
end

function get_bas_mag_pend(res)
    data, file = produce_or_load(
        datadir("basins"),
        @dict(res),
        _get_bas_mag_pend;
        prefix = "basin_mag"
    )
    @unpack basins = data
    return basins
end


# ── Newton map ────────────────────────────────────────────────────────────────

function newton_map!(dz, z, p, n)
    f(x)  = x^p[1] - 1
    df(x) = p[1] * x^(p[1]-1)
    z1    = z[1] + im*z[2]
    z1    = z1 - f(z1)/df(z1)
    dz[1] = real(z1)
    dz[2] = imag(z1)
    return
end

function get_basins_newt_map(res)
    ds     = DeterministicIteratedMap(newton_map!, [0.1, 0.2], [3])
    xg     = range(-1.0, 1.0, length=res)
    yg     = range(-1.0, 1.0, length=res)
    mapper = AttractorsViaRecurrences(ds, (xg, yg))
    basins, attractors = basins_of_attraction(mapper, (xg, yg))
    return basins
end
