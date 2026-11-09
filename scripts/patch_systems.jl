# Patch metrics of the dynamical systems, the Julia sets (rough boundary, no
# proliferation) and the carpets. Table patch_systems.tex, figure patch_systems.pdf.
include(joinpath(@__DIR__, "patch_common.jl"))

cases = [
    ("newton",     "Newton \$z^3-1\$",      get_basins_newton, [200, 300, 500, 700, 1000, 1400, 2000]),
    ("magnetic",   "Magnetic pendulum",     get_basins_magnetic, [100, 150, 200, 300, 400, 500]),
    ("feudel",     "Feudel map",            get_basins_feudel, [100, 150, 200, 300, 400, 500]),
    ("duffing",    "Duffing (Wada)",        N -> get_basins_duffing(DUFFING_F, N), [100, 150, 200, 300, 400, 500, 700, 1000]),
    [("julia_" * replace(n, " " => "_"), "Julia, " * n, N -> julia_basins(c, N), [200, 300, 400, 600, 800, 1200, 1600])
     for (n, c, _) in JULIA_CASES]...,
    ("carpet",     "Sierpi\\'nski carpet",  N -> sierpinski_carpet(round(Int, log(3, N))), 3 .^ (3:7)),
]

results = Dict(name => run_case(name, get, Ns) for (name, _, get, Ns) in cases)

rows = map(cases) do (name, label, _, _)
    r = results[name]; m = r.last
    [label, string(m.NA), fmt(r.d_b), fmt(r.α), fmt(m.twoK), fmt(m.β; d = 2), fmt(m.p; d = 2),
     fmt(m.γ; d = 2), fmt(m.τ; d = 2), fmt(r.α / (2m.β); d = 2), fmt(m.tri_per_v; d = 2),
     fmt(m.wada; d = 2), fmte(m.a_lr), fmte(m.a_tb), fmt(r.kmax_exp; d = 2)]
end
write_table("patch_systems.tex",
    "Patch metrics of the basins and reference sets. \$d_b\$ and \$\\alpha\$: exponents of the
    number of boundary pixels and of patches with the resolution; the other columns at the finest
    grid. \$a^*_\\leftrightarrow\$, \$a^*_\\updownarrow\$: bottleneck for left--right and top--bottom
    crossings (fraction of the map); \$k_{\\max}\\sim|V|^{x}\$: growth of the largest degree.",
    "tab:systems",
    ["System", "\$N_A\$", "\$d_b\$", "\$\\alpha\$", "\$2K\$", "\$\\beta\$", "\$p\$", "\$\\gamma\$",
     "\$\\tau\$", "\$\\alpha/2\\beta\$", "tri/\$|V|\$", "Wada", "\$a^*_\\leftrightarrow\$",
     "\$a^*_\\updownarrow\$", "\$x\$"], rows)

# (a) bottleneck against N, (b) degree CCDF at the finest grid with the predicted tail
fig = Figure(size = (640, 280))
ax1 = Axis(fig[1, 1]; xlabel = L"N", ylabel = L"a^*_\leftrightarrow", xscale = log10, yscale = log10,
           title = L"(a)\ \mathrm{bottleneck}")
for (name, label, _, _) in cases[1:4]
    r = results[name]
    scatterlines!(ax1, r.Ns, [max(m.a_lr, 1e-7) for m in r.per_N];
                  label = replace(label, "\$" => "", "z^3-1" => "z³−1"))
end
axislegend(ax1; position = :rc, framevisible = false, labelsize = 9)
ax2 = Axis(fig[1, 2]; xlabel = L"k", ylabel = L"P(k'\geq k)", xscale = log10, yscale = log10,
           title = L"(b)\ \mathrm{degree\ distribution}")
for (k, (name, label, get, Ns)) in enumerate(cases[1:3])
    g = PatchMap(get(last(Ns))).g; deg = sort(degree(g)); deg = deg[deg .> 0]
    ks = unique(deg); ccdf = [count(≥(x), deg) / length(deg) for x in ks]
    c = Makie.wong_colors()[k]
    lines!(ax2, ks, ccdf; color = c, label = replace(label, "\$" => "", "z^3-1" => "z³−1"))
    τp = results[name].α / (2results[name].last.β)   # predicted tail α/(2β)
    x0 = 6; y0 = count(≥(x0), deg) / length(deg)
    lines!(ax2, [x0, maximum(ks)], [y0, y0 * (maximum(ks) / x0)^(-τp)]; color = c, linestyle = :dash)
end
axislegend(ax2; position = :lb, framevisible = false, labelsize = 9)
save(joinpath(FIGDIR, "patch_systems.pdf"), fig; pt_per_unit = 1)
