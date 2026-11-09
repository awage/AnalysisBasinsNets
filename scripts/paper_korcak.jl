# Korčák analysis of the BAG paper (draft/main.tex): vertex scaling α, single-grid
# Korčák dimension 2K, boundary-pixel exponent d_b and the Wada averages. Reads the results stored by
# patch_systems.jl and patch_saupe.jl (data/patch_results); only the mean Wada
# count w̄ is recomputed, on the finest grid of each system.
# Writes draft/tables/korcak.tex and draft/figures/korcak.pdf.
include(joinpath(@__DIR__, "patch_common.jl"))

systems = [
    ("newton",   "Newton \$z^3-1\$",   get_basins_newton),
    ("magnetic", "Magnetic pendulum",  get_basins_magnetic),
    ("feudel",   "Feudel map",         get_basins_feudel),
    ("duffing",  "Duffing oscillator", N -> get_basins_duffing(DUFFING_F, N)),
]

rows = map(systems) do (name, label, get)
    r = load_patch(name); m = r.last
    wbar = wada_mean(PatchMap(get(last(r.Ns))))
    wr = extrema(s.wada for s in r.per_N)
    @printf("%-9s α=%.3f d_b=%.3f 2K=%.3f w̄=%.3f Wada∈[%.3f, %.3f]\n",
            name, r.α, r.d_b, m.twoK, wbar, wr...)
    [label, string(m.NA), fmt(r.d_b), fmt(r.α),
     fmt(m.twoK), @sprintf("%.2f--%.2f", wr...), fmt(wbar; d = 2)]
end
write_table("korcak.tex",
    "Korčák analysis of the BAG. \$d_b\$ and \$\\alpha\$: exponents of the number of boundary
    pixels and of vertices with the resolution. \$2K\$: Korčák
    dimension fitted on the finest grid. Wada: range over \$N\$ of the fraction of boundary
    components adjacent to all other basins; \$\\bar w\$: mean number of other basins adjacent to
    a boundary component, on the finest grid.",
    "tab:korcak",
    ["System", "\$N_A\$", "\$d_b\$", "\$\\alpha\$", "\$2K\$", "Wada", "\$\\bar w\$"],
    rows; wide = false)

# Saupe family: boundary-pixel exponent, vertex exponent and single-grid 2K against d_F.
hs = sort(collect(keys(SAUPE_DF)))
saupe = [load_patch("saupe_h=$(h)") for h in hs]
dF = [SAUPE_DF[h] for h in hs]

fig = Figure(size = (300, 280))
ax = Axis(fig[1, 1]; xlabel = L"d_F\ \mathrm{(Saupe)}", ylabel = "exponent", aspect = 1)
lines!(ax, [1.2, 1.72], [1.2, 1.72]; color = :gray, linestyle = :dash)
scatter!(ax, dF, [r.d_b for r in saupe]; label = L"d_b", marker = :rect)
scatter!(ax, dF, [r.α for r in saupe]; label = L"\alpha")
scatter!(ax, dF, [r.last.twoK for r in saupe]; label = L"2K", marker = :utriangle)
axislegend(ax; position = :lt, framevisible = false)
save(joinpath(FIGDIR, "korcak.pdf"), fig; pt_per_unit = 1)
println("wrote ", joinpath(FIGDIR, "korcak.pdf"))
