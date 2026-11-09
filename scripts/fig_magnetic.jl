# Figure 6: magnetic pendulum with three magnets, damping 0.05.
include(joinpath(@__DIR__, "paper_common.jl"))

res = analyze_system("magnetic", get_basins_magnetic, NS)

N = 300; g = range(MAG_BOX..., length = N)
mags = mag_magnets()
plot_basins("magnetic_basins.pdf", g, g, get_basins_magnetic(N);
            title = "Magnetic pendulum, N = $N",
            extra! = ax -> scatter!(ax, first.(mags), last.(mags); marker = :star5,
                                    color = :black, markersize = 12))
plot_scaling("magnetic_scaling.pdf", [res]; title = "Magnetic pendulum BAG vertex scaling")
