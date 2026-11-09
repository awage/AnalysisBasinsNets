# Reading and rasterising the MFE50 Mediterranean mixed-forest layer
# (../data/mediterranean_mixed_forests/TodosMed_Ecorreg.*) without GIS
# libraries: a minimal ESRI shapefile (polygon) and dBase reader, and a
# scanline rasteriser with the even–odd rule. Coordinates are UTM 30N (ED50),
# in metres.

const FOREST_DIR = joinpath(dirname(dirname(@__DIR__)), "data", "mediterranean_mixed_forests")
const FOREST_BASE = joinpath(FOREST_DIR, "TodosMed_Ecorreg")

struct ShpPolygon
    bbox::NTuple{4,Float64}          # xmin, ymin, xmax, ymax
    parts::Vector{Int32}             # 0-based start index of each ring
    x::Vector{Float64}
    y::Vector{Float64}
end

"""
    read_shp_polygons(path) -> Vector{ShpPolygon}

Polygons (shape type 5) of an ESRI shapefile, in file order (= dBase row order).
"""
function read_shp_polygons(path = FOREST_BASE * ".shp")
    data = read(path)
    rd(T, off) = ltoh(reinterpret(T, data[off+1:off+sizeof(T)])[1])
    rdbe(off) = ntoh(reinterpret(Int32, data[off+1:off+4])[1])
    rd(Int32, 32) == 5 || error("expected a polygon shapefile")
    polys = ShpPolygon[]
    off = 100
    while off < length(data)
        len = 2 * rdbe(off + 4)                 # content length in bytes
        c = off + 8
        if rd(Int32, c) == 5
            bbox = ntuple(i -> rd(Float64, c + 4 + 8(i - 1)), 4)
            np = rd(Int32, c + 36); npt = rd(Int32, c + 40)
            parts = [rd(Int32, c + 44 + 4i) for i in 0:np-1]
            p0 = c + 44 + 4np
            xy = reinterpret(Float64, data[p0+1:p0+16npt])
            push!(polys, ShpPolygon(bbox, parts, xy[1:2:end], xy[2:2:end]))
        else                                     # null shape
            push!(polys, ShpPolygon((NaN, NaN, NaN, NaN), Int32[], Float64[], Float64[]))
        end
        off = c + len
    end
    return polys
end

"""
    read_dbf_columns(path, names) -> Dict{String,Vector{String}}

Selected columns of a dBase III file as stripped strings.
"""
function read_dbf_columns(path, names)
    io = open(path)
    h = read(io, 32)
    n = Int(reinterpret(UInt32, h[5:8])[1]); hl = Int(reinterpret(UInt16, h[9:10])[1])
    rl = Int(reinterpret(UInt16, h[11:12])[1])
    fields = Tuple{String,Int,Int}[]; pos = 1     # (name, offset in record, length)
    while true
        d = read(io, 32); d[1] == 0x0d && break
        z = findfirst(==(0x00), d[1:11]); name = String(d[1:(z === nothing ? 11 : z - 1)])
        push!(fields, (name, pos, Int(d[17]))); pos += Int(d[17])
    end
    seek(io, hl); raw = read(io, n * rl); close(io)
    out = Dict{String,Vector{String}}()
    for nm in names
        f = findfirst(t -> t[1] == nm, fields); f === nothing && error("no field $nm")
        _, o, L = fields[f]
        out[nm] = [strip(latin1(raw[(i-1)*rl+o+1:(i-1)*rl+o+L])) for i in 1:n]
    end
    return out
end

# dBase strings are Latin-1: one byte per character.
latin1(v::AbstractVector{UInt8}) = String(map(Char, v))

"""
    rasterize(polys, labels, cell; bbox) -> Matrix{Int32}

Raster with square cells of side `cell` (metres) over `bbox` = (xmin, ymin,
xmax, ymax). A cell takes the label of the polygon that contains its centre
(even–odd rule over all rings); 0 elsewhere. Later polygons overwrite earlier
ones. Rows are x (west→east), columns y (south→north).
"""
function rasterize(polys, labels, cell; bbox)
    xmin, ymin, xmax, ymax = bbox
    nx = ceil(Int, (xmax - xmin) / cell); ny = ceil(Int, (ymax - ymin) / cell)
    R = zeros(Int32, nx, ny)
    for (p, lab) in zip(polys, labels)
        (lab == 0 || isempty(p.x)) && continue
        # columns whose centre y lies inside the polygon's bounding box
        j0 = max(1, ceil(Int, (p.bbox[2] - ymin) / cell + 0.5))
        j1 = min(ny, floor(Int, (p.bbox[4] - ymin) / cell + 0.5))
        j0 > j1 && continue
        crossings = [Float64[] for _ in j0:j1]
        bounds = vcat(p.parts, Int32(length(p.x)))
        for r in 1:length(p.parts)
            a, b = bounds[r] + 1, bounds[r+1]
            for k in a:b-1
                x1, y1, x2, y2 = p.x[k], p.y[k], p.x[k+1], p.y[k+1]
                y1 == y2 && continue
                lo, hi = minmax(y1, y2)
                ja = max(j0, ceil(Int, (lo - ymin) / cell + 0.5))
                jb = min(j1, ceil(Int, (hi - ymin) / cell + 0.5) - 1)   # half-open [lo, hi)
                for j in ja:jb
                    yc = ymin + (j - 0.5) * cell
                    push!(crossings[j-j0+1], x1 + (yc - y1) * (x2 - x1) / (y2 - y1))
                end
            end
        end
        for (jj, xs) in enumerate(crossings)
            length(xs) < 2 && continue
            sort!(xs); j = j0 + jj - 1
            for m in 1:2:length(xs)-1
                ia = max(1, ceil(Int, (xs[m] - xmin) / cell + 0.5))
                ib = min(nx, floor(Int, (xs[m+1] - xmin) / cell + 0.5))
                for i in ia:ib
                    R[i, j] = lab
                end
            end
        end
    end
    return R
end

"""
    forest_layer(field) -> (polys, labels, names, area_ha)

Polygons of the layer with an integer label per distinct value of the dBase
column `field` (1, 2, … in order of decreasing frequency) and the polygon areas.
"""
function forest_layer(field)
    polys = read_shp_polygons()
    cols = read_dbf_columns(FOREST_BASE * ".dbf", [field, "Area_ha"])
    vals = cols[field]
    counts = Dict{String,Int}(); for v in vals; counts[v] = get(counts, v, 0) + 1; end
    names = sort(collect(keys(counts)), by = v -> -counts[v])
    id = Dict(v => i for (i, v) in enumerate(names))
    labels = Int32[id[v] for v in vals]
    return polys, labels, names, parse.(Float64, cols["Area_ha"])
end

layer_bbox(polys) = (minimum(p.bbox[1] for p in polys if !isempty(p.x)),
                     minimum(p.bbox[2] for p in polys if !isempty(p.x)),
                     maximum(p.bbox[3] for p in polys if !isempty(p.x)),
                     maximum(p.bbox[4] for p in polys if !isempty(p.x)))
