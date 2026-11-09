# Relaxed Newton family of Saupe (1987): BAG metrics against the published
# box-counting dimension d_F(h). Table patch_saupe.tex, figure patch_saupe.pdf.
include(joinpath(@__DIR__, "patch_common.jl"))

const SAUPE_NS = [250, 354, 500, 707, 1000, 1414, 2000]
hs = sort(collect(keys(SAUPE_DF)))

results = map(hs) do h
    run_case("saupe_h=$(h)", N -> relaxed_newton_basins(h, N), SAUPE_NS)
end

rows = map(zip(hs, results)) do (h, r)
    m = r.last
    [fmt(h; d = 1), fmt(SAUPE_DF[h]), fmt(r.d_b), fmt(r.α), fmt(m.twoK), fmt(m.β; d = 2),
     fmt(m.p; d = 2), fmt(m.γ; d = 2), fmt(m.τ; d = 2), fmt(r.α / (2m.β); d = 2),
     fmt(m.tri_per_v; d = 2), fmte(m.a_lr)]
end
write_table("patch_saupe.tex",
    "Relaxed Newton family \$N_h\$ of Saupe~\\cite{saupe1987}. \$d_F\$: box-counting dimension of the
    Julia set (Saupe, Table IV); \$d_b\$: exponent of the number of boundary pixels; \$\\alpha\$:
    vertex scaling, both for \$N=250\\ldots2000\$. The remaining columns are computed on the
    \$N=2000\$ grid: Korčák dimension \$2K\$, contact exponents \$\\beta=p\\gamma\$, degree tail
    \$\\tau\$ and its prediction \$\\alpha/2\\beta\$, triangles per patch and bottleneck \$a^*\$.",
    "tab:saupe",
    ["\$h\$", "\$d_F\$", "\$d_b\$", "\$\\alpha\$", "\$2K\$", "\$\\beta\$", "\$p\$", "\$\\gamma\$",
     "\$\\tau\$", "\$\\alpha/2\\beta\$", "tri/\$|V|\$", "\$a^*\$"], rows)

fig = Figure(size = (620, 270))
ax1 = Axis(fig[1, 1]; xlabel = L"d_F\ \mathrm{(Saupe)}", ylabel = "exponent",
           title = "(a) dimensions", aspect = 1)
dF = [SAUPE_DF[h] for h in hs]
lines!(ax1, [1.2, 1.72], [1.2, 1.72]; color = :gray, linestyle = :dash)
scatter!(ax1, dF, [r.d_b for r in results]; label = L"d_b", marker = :rect)
scatter!(ax1, dF, [r.α for r in results]; label = L"\alpha")
scatter!(ax1, dF, [r.last.twoK for r in results]; label = L"2K", marker = :utriangle)
axislegend(ax1; position = :lt, framevisible = false)
ax2 = Axis(fig[1, 2]; xlabel = L"h", ylabel = "exponent", title = "(b) contact exponents")
for (sym, lab, mk) in ((:β, L"\beta", :circle), (:p, L"p", :rect), (:γ, L"\gamma", :utriangle))
    scatterlines!(ax2, hs, [getfield(r.last, sym) for r in results]; label = lab, marker = mk)
end
axislegend(ax2; position = :lt, framevisible = false)
save(joinpath(FIGDIR, "patch_saupe.pdf"), fig; pt_per_unit = 1)
