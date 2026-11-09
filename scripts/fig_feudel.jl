# Figure 7: quasiperiodically forced map of Feudel et al. (1998).
include(joinpath(@__DIR__, "paper_common.jl"))

att = feudel_attractors()
println("Feudel attractors: ", length(att), "  sizes: ", [length(a) for a in values(att)])

res = analyze_system("feudel", get_basins_feudel, NS)

N = 200; g = range(0, 1, length = N)
plot_basins("feudel_basins.pdf", g, g, get_basins_feudel(N);
            xlabel = L"u", ylabel = L"x",
            title = "Feudel map (r = $(FEUDEL_PARAMS.r), a = $(FEUDEL_PARAMS.a))")
plot_scaling("feudel_scaling.pdf", [res]; title = "Feudel BAG vertex scaling")
