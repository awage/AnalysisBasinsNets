# Figure 3: Newton map for z^3 - 1.
include(joinpath(@__DIR__, "paper_common.jl"))

res = analyze_system("newton", get_basins_newton, vcat(NS_LARGE, [1400, 2000]))

# Grid-free reference: uncertainty exponent from random pairs of initial
# conditions (McDonald et al. 1985), d = 2 - α_u.
α_u, = grebogi_exponent((x, y) -> newton_root(complex(x, y)), NEWTON_BOX)
res["alpha_u_grebogi"] = α_u
res["d_grebogi"] = 2 - α_u
@printf("newton  grid-free uncertainty exponent α_u = %.3f  →  d = %.3f\n", α_u, 2 - α_u)
wsave(datadir("paper_results", "newton.jld2"), res)

N = 200; r = range(NEWTON_BOX..., length = N)
plot_basins("newton_basins.pdf", r, r, get_basins_newton(N);
            xlabel = L"\mathrm{Re}\, z", ylabel = L"\mathrm{Im}\, z",
            title = "Newton basins (z³ − 1), N = $N")
plot_scaling("newton_scaling.pdf", [res]; title = "Newton BAG vertex scaling")
