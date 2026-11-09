# Patch metrics of the MFE50 Mediterranean mixed forests (Spain): the vector
# layer is rasterised at cell sizes from 100 m to 1.6 km, the non-forest matrix
# (label 0) is removed from the graph, and the metrics of draft/patchiness_note.tex
# are computed for two classifications: TIPO_BOSQU (3 classes: broadleaf,
# conifer, mixed) and Bosqumix (24 mixed-forest types).
# Results: data/patch_results/forest_<field>.jld2, table draft/tables/patch_forest.tex.
include(joinpath(@__DIR__, "patch_common.jl"))
include(srcdir("forest_map.jl"))

const FOREST_CELLS = [100.0, 141.0, 200.0, 283.0, 400.0, 566.0, 800.0, 1131.0, 1600.0]

function forest_sweep(field)
    polys, labels, names, _ = forest_layer(field)
    bb = layer_bbox(polys)
    ms = map(FOREST_CELLS) do cell
        t = @elapsed begin
            R = rasterize(polys, labels, cell; bbox = bb)
            m = fixed_scale_metrics(R; background = 0)
        end
        @printf("%-10s cell %5.0f m |V|=%7d 2K=%.3f β=%.3f p=%.3f γ=%.3f τ=%.2f tri/V=%.3f wada=%.3f kmax=%d hops=%.1f (%.0fs)\n",
                field, cell, m.nv, m.twoK, m.β, m.p, m.γ, m.τ, m.tri_per_v, m.wada, m.kmax, m.hops, t)
        flush(stdout)
        m
    end
    nv = [m.nv for m in ms]
    local_α = -diff(log.(nv)) ./ diff(log.(FOREST_CELLS))
    α = -ols_slope(FOREST_CELLS, nv)
    d_b = -ols_slope(FOREST_CELLS, [m.nb for m in ms])
    @printf("%-10s α = %.3f  d_b = %.3f  local α (fine→coarse) = %s\n", field, α, d_b,
            string(round.(local_α, digits = 2)))
    res = (; field, names, cells = FOREST_CELLS, per_cell = ms, α, d_b, local_α)
    save_patch("forest_$(field)", res)
    return res
end

results = [forest_sweep(f) for f in ("TIPO_BOSQU", "Bosqumix")]

rows = map(results) do r
    m = r.per_cell[1]        # finest cell
    [r.field == "TIPO_BOSQU" ? "Forest type" : "Mixed-forest type", string(m.NA), string(m.nv),
     fmt(100m.isolated; d = 0), fmt(m.mean_k; d = 2), fmt(r.local_α[1]; d = 2),
     fmt(r.local_α[end]; d = 2), fmt(r.d_b; d = 2), fmt(m.β; d = 2), fmt(m.p; d = 2),
     fmt(m.tri_per_v; d = 3), fmt(m.wada; d = 2)]
end
write_table("patch_forest.tex",
    "MFE50 Mediterranean mixed forests of Spain rasterised at 100\\,m--1.6\\,km, non-forest matrix
    removed from the graph. Patches, isolated patches (no contact with another forest class), mean
    degree, \$\\beta\$, \$p\$, triangles per patch and Wada fraction at 100\\,m. The vertex exponent
    is not constant: local \$\\alpha\$ between 100 and 141\\,m and between 1.1 and 1.6\\,km.
    \$d_b\$: boundary-pixel exponent over the whole range.",
    "tab:forest",
    ["Classes", "\$N_A\$", "\$|V|\$", "isolated (\\%)", "\$\\langle k\\rangle\$",
     "\$\\alpha_{\\rm fine}\$", "\$\\alpha_{\\rm coarse}\$", "\$d_b\$", "\$\\beta\$", "\$p\$",
     "tri/\$|V|\$", "Wada"], rows)

# ── Korčák profile of the patch areas and comparison with the local vertex exponent ──
# Patch areas (ha) of the forest-type map at 100 m (1 pixel = 1 ha). The local
# Korčák dimension 2K(a) is the slope of N(A ≥ a) over one decade of area around
# a; the resolution sweep should give α(c) ≈ 2K(a = c²).
polys, labels, _, _ = forest_layer("TIPO_BOSQU")
pm = PatchMap(rasterize(polys, labels, 100.0; bbox = layer_bbox(polys)); background = 0)
A = Float64.(pm.area)
us = sort(unique(A)); Nge = [count(≥(u), A) for u in us]
function local_2K(a; halfwidth = 0.5)          # window of ± halfwidth decades
    s = (us .≥ a / 10^halfwidth) .& (us .≤ a * 10^halfwidth)
    count(s) < 3 && return NaN
    x = log.(us[s]); y = log.(Nge[s]); -2cov(x, y) / var(x)
end
decades = [(4, 40), (40, 400), (400, 4000), (4000, 40000)]
for (a0, a1) in decades
    @printf("areas %5d–%5d ha: 2K = %.2f (%d patches ≥ %d ha)\n", a0, a1, local_2K(sqrt(a0 * a1)),
            count(≥(a0), A), a0)
end
save_patch("forest_korcak_profile", (; us, Nge, decades,
           twoK = [local_2K(sqrt(a0 * a1)) for (a0, a1) in decades]))

fig = Figure(size = (640, 280))
ax1 = Axis(fig[1, 1]; xlabel = L"a\ \mathrm{(ha)}", ylabel = L"N(A\geq a)", xscale = log10,
           yscale = log10, title = L"(a)\ \mathrm{Kor\check{c}\acute{a}k\ plot,\ 100\ m}")
lines!(ax1, us, Nge; color = :black)
for (sl, c) in ((0.5, :gray), (1.0, :gray40))   # reference slopes K = 0.5 and 1
    lines!(ax1, [10.0, 1000.0], 2e4 .* ([10.0, 1000.0] ./ 10) .^ (-sl); color = c, linestyle = :dash,
           label = "K = $sl")
end
axislegend(ax1; position = :lb, framevisible = false)
ax2 = Axis(fig[1, 2]; xlabel = L"\mathrm{cell}\ c\ \mathrm{(m)}", ylabel = "local exponent",
           xscale = log10, title = L"(b)\ \alpha(c)\ \mathrm{and}\ 2K(a=c^2)")
for (r, mk) in zip(results, (:circle, :rect))
    cm = sqrt.(r.cells[1:end-1] .* r.cells[2:end])
    scatterlines!(ax2, cm, r.local_α; marker = mk,
                  label = r.field == "TIPO_BOSQU" ? L"\alpha,\ 3\ \mathrm{classes}" : L"\alpha,\ 24\ \mathrm{classes}")
end
cs = 10 .^ range(log10(150), log10(3000), length = 25)
lines!(ax2, cs, [local_2K(c^2 / 1e4) for c in cs]; color = :black, linestyle = :dash, label = L"2K(c^2)")
axislegend(ax2; position = :lt, framevisible = false)
save(joinpath(FIGDIR, "patch_forest.pdf"), fig; pt_per_unit = 1)
