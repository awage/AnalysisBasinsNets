# Neutral random quadtree mosaics: calibration of α, β, τ and the bottleneck a*.
# Table patch_neutral.tex, figure patch_neutral.pdf.
include(joinpath(@__DIR__, "patch_common.jl"))

qs = [0.35, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95]
ds = 7:11
modes = [(0, "unique"), (3, "3 colours")]

results = Dict()
for q in qs, (cl, mode) in modes
    results[(q, cl)] = run_case("quadtree_q=$(q)_colors=$(cl)",
                                N -> quadtree_mosaic(q, round(Int, log2(N)); colors = cl), 2 .^ ds)
end

βtheory(q) = max(0, log2(2q) / 2)
rows = map(Iterators.product(qs, modes) |> collect |> vec) do (q, (cl, mode))
    r = results[(q, cl)]; m = r.last
    [mode, fmt(q; d = 2), fmt(r.α), fmt(log2(4q)), fmt(m.β; d = 2), fmt(βtheory(q); d = 2),
     fmt(m.p; d = 2), fmt(m.γ; d = 2), fmt(m.τ; d = 2), fmt(r.α / (2m.β); d = 2),
     fmt(m.tri_per_v; d = 2), fmte(m.a_lr)]
end
write_table("patch_neutral.tex",
    "Neutral random quadtree mosaics (a cell splits into four with probability \$q\$), rasterised
    at \$N=2^7\\ldots2^{11}\$. With unique labels every leaf is a patch; with three random colours
    equal neighbours merge. Expected \$\\alpha=\\log_2 4q\$ and, for unique labels,
    \$\\beta=\\max(0,\\tfrac12\\log_2 2q)\$. Fixed-scale columns at \$N=2048\$.",
    "tab:neutral",
    ["labels", "\$q\$", "\$\\alpha\$", "\$\\log_2 4q\$", "\$\\beta\$", "theory", "\$p\$",
     "\$\\gamma\$", "\$\\tau\$", "\$\\alpha/2\\beta\$", "tri/\$|V|\$", "\$a^*\$"], rows)

fig = Figure(size = (620, 270))
ax1 = Axis(fig[1, 1]; xlabel = L"q", ylabel = L"\beta", title = L"(a)\ \beta\ \mathrm{vs\ split\ probability}")
lines!(ax1, 0.3:0.01:1.0, βtheory.(0.3:0.01:1.0); color = :gray, linestyle = :dash, label = "theory (unique)")
for (cl, mode) in modes
    scatterlines!(ax1, qs, [results[(q, cl)].last.β for q in qs]; label = mode)
end
axislegend(ax1; position = :lt, framevisible = false)
ax2 = Axis(fig[1, 2]; xlabel = L"N", ylabel = L"a^*", xscale = log10, yscale = log10,
           title = L"(b)\ \mathrm{bottleneck,\ 3\ colours}")
for q in qs
    r = results[(q, 3)]
    scatterlines!(ax2, r.Ns, [max(m.a_lr, 1e-7) for m in r.per_N]; label = "q = $q")
end
axislegend(ax2; position = :lb, framevisible = false, labelsize = 9)
save(joinpath(FIGDIR, "patch_neutral.pdf"), fig; pt_per_unit = 1)
