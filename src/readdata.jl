using MAT
using VideoIO
using Images

#const datadir = "D:/"
const datadir = normpath("$(@__DIR__)/../data")

function prefix( cameradir, datadir = datadir )
    files = readdir(joinpath( datadir, cameradir ))
    path, dir = splitdir(cameradir)
    matches = [match( Regex(".*($(dir).*?)\\d+.mat", "i"), f ) for f in files]
    prefixes = Set( m[1] for m in matches if m !== nothing )
    length(prefixes) > 1 && error("More than one prefix in dir $(dir)")
    first(prefixes)
end

function filepath_f( cameradir, datadir = datadir )
    pr = prefix( cameradir, datadir )
    (type, idx) -> "$datadir/$cameradir/$type$pr$(string(idx;pad=4)).mat"
end

#filetypes = [(:coords, "corrd", "x_y_coor"),
filetypes = [(:coords, "coord", "x_y_coor"),
             (:cropped, "short", "mydata"),
             #(:trace, "Trajectory", "Trace"),
             (:size, "WormSize", "big")]

for (sym, prefix, varname) in filetypes
    @eval function $(Symbol(:read_,sym))(path_f, idx)
        tr = matopen(path_f($prefix,idx)) do file
            read(file,$varname)
        end
    end
end

# read old-format cropped frames
function readframes( path_f, idx )
    frames = vec( read_cropped( path_f, idx ) )
    [colorview(RGB, normedview(permuteddimsview(fr,(3,1,2)))) for fr in frames]
end

function _readall!(f, firstframe)
    out = [firstframe]
    while !eof(f)
        push!(out, read(f))
    end
    out
end
_readall!(f) = _readall!(f,read(f))

function video_prefix( cameradir, datadir = datadir )
    files = readdir(joinpath( datadir, cameradir ))
    path, dir = splitdir(cameradir)
    dir = chop(dir,tail=2)
    matches = [match( Regex("shape($(dir).*?)\\d+.mp4", "i"), f ) for f in files]
    prefixes = Set( m[1] for m in matches if m !== nothing )
    length(prefixes) > 1 && error("More than one prefix in dir $(dir)")
    first(prefixes)
end

function videopath_f( cameradir, datadir = datadir )
    pr = video_prefix( cameradir, datadir )
    idx -> "$datadir/$cameradir/shape$pr$(string(idx;pad=4)).mp4"
end

function read_video(path_f, idx)
    vid = openvideo(path_f(idx))
    try
        _readall!(vid)
    finally
        close(vid)
    end
end

nfiles(path_f) = first( i for i in Iterators.countfrom(0)
                             if !isfile(path_f("coord",i)) )
#                             if !isfile(path_f("corrd",i)) )

using DataFrames

function import_coords( cameradir, datadir = datadir; with_size=false )
    path_f = filepath_f(cameradir, datadir)
    x = Union{Float64,Missing}[]
    y = Union{Float64,Missing}[]
    with_size && (sz = Union{Float64,Missing}[])
    fileno = Int[]
    batch_boundaries = Int[]
    n = 0
    #@progress "Importing"
    for i in 0:nfiles(path_f)-1
        coords = read_coords(path_f,i)
        if size(coords,1) == 2 && size(coords,2) > 2
            coords = coords'
        end
        append!(x, @view coords[:,1])
        append!(y, @view coords[:,2])
        with_size && append!(sz, read_size(path_f,i))
        prev_n, n = n, length(x)
        resize!(fileno, n)
        fileno[prev_n+1:n] .= i
        push!(batch_boundaries, n)
    end
    # The old scripts use 1.0 when worm tracking fails. New scripts use NaN.
    mark_missings(a) = a[isequal.(a,1) .| isnan.(a)] .= missing
    mark_missings(x)
    mark_missings(y)
    with_size && mark_missings(sz)
    
    df = DataFrame( fileno = fileno, x = x, y = y )
    with_size && (df.size = sz)
    df, batch_boundaries
end

_video_index(boundaries, frameno) = 0 < frameno <= last(boundaries) ?
                                    searchsortedfirst(boundaries,frameno) :
                                    nothing

const frames_per_s = 3

# TODO pixel size per experiment? (pixels were 10μm in older experiments)
# const pixel_size_μm = 10
