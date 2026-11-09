# Notes only (draft/contact_exponents_notes.tex, Sec. "Degree tail and the Korcak law"), not used in the paper.
# Test: degree tail <-> Korcak via the edge size-ratio distribution.
# Prediction: beta = K - theta, tau = K/(K - theta), with P(A_big/A_small > r) ~ r^-theta over edges.
include(joinpath(@__DIR__, "patch_common.jl"))

hill(x, xmin) = (t = x[x .≥ xmin]; length(t) / sum(log.(t ./ xmin)))

function analyse(name, B)
    pm = PatchMap(B); g = pm.g; A = pm.area
    k = degree(g)
    K, _ = korcak_exponent(pm)
    β, p, γ, _ = contact_exponents(pm)
    τ = degree_tail(pm)
    # edge ratios and up/down degree
    r = Float64[]; up = zeros(Int, nv(g)); dn = zeros(Int, nv(g))
    for e in edges(g)
        i, j = src(e), dst(e)
        a, b = minmax(A[i], A[j])
        push!(r, b / a)
        if A[i] < A[j]; up[i] += 1; dn[j] += 1
        elseif A[j] < A[i]; up[j] += 1; dn[i] += 1 end
    end
    θ = hill(r, 10.0)
    # D(>a): edges whose larger endpoint exceeds a  ~ a^(-θ')
    bigend = [max(A[src(e)], A[dst(e)]) for e in edges(g)]
    us = exp.(range(log(10), log(1e-3 * pm.npix), length = 15))
    D = [count(>(u), bigend) for u in us]
    θD = -ols_slope(us, Float64.(D))
    s = k .> 0
    @printf("%-14s |V|=%7d <k>=%.2f  K=%.3f β=%.3f τ=%.2f | θ_hill=%.3f θ_D=%.3f | K-θ=%.3f (β)  K/(K-θ)=%.2f (τ) | <up>=%.2f  P(up≤3)=%.2f  med k=%d\n",
            name, count(s), mean(k[s]), K, β, τ, θ, θD, K - θD, K / (K - θD), mean(up[s]), mean(up[s] .≤ 3), median(k[s]))
    flush(stdout)
end

analyse("newton", get_basins_newton(2000))
analyse("magnetic", get_basins_magnetic(500))
analyse("feudel", get_basins_feudel(500))
for q in (0.5, 0.7, 0.9)
    analyse("qt3 q=$q", quadtree_mosaic(q, 11; colors = 3))
    analyse("qtU q=$q", quadtree_mosaic(q, 11; colors = 0))
end
