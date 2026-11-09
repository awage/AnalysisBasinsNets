# Reproduce the tables and figures of draft/patchiness_note.tex.
# Run with threads:  julia -t auto --project=. scripts/make_patchiness.jl
for s in ("patch_lacunarity.jl", "patch_neutral.jl", "patch_systems.jl", "patch_saupe.jl",
          "forest_patches.jl")
    @info "Running $s"
    path = joinpath(@__DIR__, s)
    Core.eval(Main, :(module $(Symbol(splitext(s)[1])) include($path) end))
end
