# Figures 1 and 2: Sierpiński triangle and carpet as artificial basins.
include(joinpath(@__DIR__, "paper_common.jl"))

ms_tri = 2:10
ms_car = 2:7
tri = analyze_system("sierpinski_triangle", N -> sierpinski_triangle(Int(log2(N))),
                     [2^m for m in ms_tri]; fit_min = 0)
car = analyze_system("sierpinski_carpet", N -> sierpinski_carpet(round(Int, log(3, N))),
                     [3^m for m in ms_car]; fit_min = 0)

# Exact vertex count for the triangle construction, to compare with Eq. (1)
for (m, s) in zip(ms_tri, tri["stats"])
    @printf("triangle m=%2d  |V|=%8d   (3^(m+1)-1)/2=%8d\n", m, s.nv, (3^(m+1) - 1) ÷ 2)
end

m = 5; B = sierpinski_triangle(m)[:, end:-1:1]; r = range(0, 1, length = 2^m)
plot_basins("sierpinski_triangle.pdf", r, r, B; xlabel = "", ylabel = "",
            title = "Sierpiński triangle (level m=$m)",
            colors = Dict(1 => colorant"#d62728", 2 => colorant"#2ca02c",
                          3 => colorant"#ff7f0e", 4 => colorant"#1f4fd6"))

m = 4; B = sierpinski_carpet(m); r = range(0, 1, length = 3^m)
plot_basins("sierpinski_carpet.pdf", r, r, B; xlabel = "", ylabel = "",
            title = "Sierpiński carpet (level m=$m)",
            colors = Dict(1 => colorant"black", 2 => colorant"white"))

tri["label"] = "triangle"; car["label"] = "carpet"
plot_scaling("sierpinski_scaling.pdf", [tri, car]; title = "Sierpiński BAG vertex scaling")
