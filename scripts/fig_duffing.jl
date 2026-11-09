# Figure 4: Duffing oscillator x'' + δx' - x + x^3 = F sin(ωt) with the Wada
# parameters of Aguirre & Sanjuán (2002): δ = 0.15, F = 0.245, ω = 1.
include(joinpath(@__DIR__, "paper_common.jl"))

res = analyze_system("duffing", N -> get_basins_duffing(DUFFING_F, N), NS_LARGE)

N = 300; g = range(DUFFING_BOX..., length = N)
att = duffing_attractors(DUFFING_F)
plot_basins("duffing_basins.pdf", g, g, get_basins_duffing(DUFFING_F, N);
            xlabel = L"x", ylabel = L"\dot{x}",
            title = "Duffing, F = $(DUFFING_F), N = $N",
            extra! = ax -> for a in values(att)
                scatter!(ax, a[:, 1], a[:, 2]; color = :black, markersize = 5)
            end)
plot_scaling("duffing_scaling.pdf", [res]; title = "Duffing BAG vertex scaling")
