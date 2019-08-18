using MAT
using VideoIO

#const datadir = "D:/"
const datadir = normpath("$(@__DIR__)/../data")

function prefix( cameradir, datadir = datadir )
    files = readdir(joinpath( datadir, cameradir ))
    path, dir = splitdir(cameradir)
    matches = [match( Regex(".*($(dir).*?)\\d+.mat"), f ) for f in files]
    prefixes = Set( m[1] for m in matches if m !== nothing )
    length(prefixes) > 1 && error("More than one prefix in dir $(dir)")
    first(prefixes)
end

function filepath_f( cameradir, datadir = datadir )
    pr = prefix( cameradir, datadir )
    (type, idx) -> "$datadir/$cameradir/$type$pr$(string(idx;pad=4)).mat"
end

filetypes = [(:coords, "corrd", "x_y_coor"),
             (:cropped, "short", "mydata"),
             (:trace, "Trajectory", "Trace"),
             (:size, "WormSize", "big")]

for (sym, prefix, varname) in filetypes
    @eval function $(Symbol(:read_,sym))(path_f, idx)
        tr = matopen(path_f($prefix,idx)) do file
            read(file,$varname)
        end
    end
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
    matches = [match( Regex(".*($(dir).*?)\\d+.mp4"), f ) for f in files]
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
                             if !isfile(path_f("corrd",i)) )

using DataFrames

function import_coords( cameradir, datadir = datadir )
    path_f = filepath_f(cameradir, datadir)
    x = Union{Float64,Missing}[]
    y = Union{Float64,Missing}[]
    fileno = Int[]
    batch_boundaries = Int[]
    n = 0
    #@progress "Importing"
    for i in 0:nfiles(path_f)-1
        array = read_coords(path_f,i)
        append!(x, @view array[:,1])
        append!(y, @view array[:,2])
        prev_n, n = n, length(x)
        resize!(fileno, n)
        fileno[prev_n+1:n] .= i
        push!(batch_boundaries, n)
    end
    # (1,1) is used when worm tracking fails
    x[isequal.(x,1)] .= missing
    y[isequal.(y,1)] .= missing
    df = DataFrame( fileno = fileno, x = x, y = y )
    df, batch_boundaries
end

const frames_per_s = 3

# TODO pixel size per experiment? (pixels were 10μm in older experiments)
# const pixel_size_μm = 10
