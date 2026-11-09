# Two carpets with the same dimension D = log 40 / log 7 and different gap
# arrangement: α is the same, the prefactor of |V| and the lacunarity differ.
# Table patch_lacunarity.tex.
include(joinpath(@__DIR__, "patch_common.jl"))

D = log(40) / log(7)
rows = map((:block, :dust)) do kind
    r = run_case("carpet7_$(kind)", N -> carpet7(kind, round(Int, log(7, N))), 7 .^ (1:4))
    nvs = [m.nv for m in r.per_N]
    local_α = diff(log.(nvs)) ./ diff(log.(r.Ns))
    B = carpet7(kind, 4)
    Λ = [gliding_box_lacunarity(B .== 2, s) for s in (7, 49, 343)]
    holes = PatchMap(B).area; holes = sort(holes)[1:end-1]   # drop the filled set
    @printf("%-6s local α = %s  |V|/N^D = %.3f  Λ = %s\n", kind, string(round.(local_α, digits = 3)),
            nvs[end] / r.Ns[end]^D, string(round.(Λ, digits = 3)))
    [string(kind), fmt(local_α[end]), fmt(D), fmt(nvs[end] / r.Ns[end]^D), fmt(mean(holes); d = 1),
     fmt(Λ[1]), fmt(Λ[2]), fmt(Λ[3])]
end
write_table("patch_lacunarity.tex",
    "Carpets on a \$7\\times7\$ generator keeping 40 cells (\$D=\\log 40/\\log 7\$): \\emph{block}
    removes the central \$3\\times3\$ square, \\emph{dust} nine isolated cells. Local \$\\alpha\$
    between levels 3 and 4, prefactor \$|V|/N^D\$, mean hole area (pixels, level 4) and gliding-box
    lacunarity \$\\Lambda(r)\$ of the filled set.",
    "tab:lacunarity",
    ["carpet", "\$\\alpha\$", "\$D\$", "\$|V|/N^D\$", "\$\\langle A_{\\rm hole}\\rangle\$",
     "\$\\Lambda(7)\$", "\$\\Lambda(49)\$", "\$\\Lambda(343)\$"], rows; wide = false)
