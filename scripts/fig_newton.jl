# Fig. 2 of the paper: Newton's method for z^3 - 1 on [-1, 1]^2.
# (a) basins at N = 200, (b) vertex scaling fitted over N ≥ 100.
# Also prints the grid-free uncertainty exponent quoted in the text.
include(joinpath(@__DIR__, "paper_common.jl"))

const NEWTON_NS = [50, 75, 100, 150, 200, 300, 400, 500, 700, 1000, 1400, 2000]
const NFIT_MIN = 100   # below, |V| is a few tens of vertices

res = analyze_system("newton", get_basins_newton, NEWTON_NS; merge_test = true)
e = exponents(res; Nmin = NFIT_MIN)
@printf("newton  α = %.3f (R² = %.4f, N ≥ %d)\n", e.α, e.R2, NFIT_MIN)

# Grid-free reference: uncertainty exponent from random pairs of initial
# conditions (McDonald et al. 1985), d = 2 - α_u.
α_u, = grebogi_exponent((x, y) -> newton_root(complex(x, y)), NEWTON_BOX)
@printf("newton  uncertainty exponent α_u = %.3f  →  d = %.3f\n", α_u, 2 - α_u)

# (a) basins
N = 200; r = range(NEWTON_BOX..., length = N)
fig = Figure(size = (300, 280))
ax = Axis(fig[1, 1]; xlabel = L"\mathrm{Re}\, z", ylabel = L"\mathrm{Im}\, z",
          title = "Newton basins (z³ − 1), N = $N", aspect = DataAspect())
image!(ax, first(r)..last(r), first(r)..last(r), basin_image(get_basins_newton(N));
       interpolate = false)
limits!(ax, first(r), last(r), first(r), last(r))
save(joinpath(FIGDIR, "newton_basins.pdf"), fig; pt_per_unit = 1)

# (b) vertex scaling; open symbols are excluded from the fit
fig = Figure(size = (320, 280))
ax = Axis(fig[1, 1]; xlabel = L"\log N", ylabel = L"\log |V|", title = "Newton BAG vertex scaling")
x = log.(res["Ns"]); y = log.([s.nv for s in res["stats"]])
out = setdiff(eachindex(x), e.fit); c = Makie.wong_colors()[1]
label = @sprintf("α = %.3f", e.α)
scatter!(ax, x[e.fit], y[e.fit]; color = c, label)
scatter!(ax, x[out], y[out]; color = :white, strokecolor = c, strokewidth = 1)
lines!(ax, x[e.fit], e.b .+ e.α .* x[e.fit]; color = c, label)
axislegend(ax; position = :lt, framevisible = false, merge = true)
save(joinpath(FIGDIR, "newton_scaling.pdf"), fig; pt_per_unit = 1)
println("wrote newton_basins.pdf, newton_scaling.pdf")
