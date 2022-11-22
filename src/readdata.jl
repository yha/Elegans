using MAT
using VideoIO
using Images
using ProgressLogging

#const datadir = "D:/"
#const datadir = normpath("$(@__DIR__)/../data")

# TODO 
# remove older "filepath_f" interface and cleanup

struct WellID
    experiment::String
    well::String
end

struct Well
    # separate `root`, `experiment` and `well` components are used in addition to full path,
    # for easier sanity checking e.g. in _import_coords_jld2
    root::String
    experiment::String
    well::String
    path::String # == joinpath(root, experiment, well)
    mat_prefix::String
    video_prefix::Union{String, Nothing}
    function Well( root, experiment, well ) 
        fullpath = joinpath( root, experiment, well )
        new(root, experiment, well, fullpath, prefix(fullpath), video_prefix(fullpath))
    end     
end

load_coords_and_size(well::Well; err_on_wrong_well=true) = load_coords_and_size( well.experiment, well.well, well.root; err_on_wrong_well )
load_video( well::Well, idx ) = if isnothing(well.video_prefix)
    @info "No videos found. Looking for short*.mat files"
    readframes( (type,i) -> _filepath( well.path, well.mat_prefix, :short, i ),
                idx )
else
    VideoIO.load( _videopath( well.path, well.video_prefix, idx ) )
end

# For sorting lists of wells
Base.isless(w1::Well, w2::Well) = let f(w) = (w.experiment, w.well, w.root)
    isless(f(w1), f(w2))
end

wellname(well::Well) = "$(well.experiment)-$(well.well)"

function prefix(welldir)
    #isnothing(datadir) || (welldir = joinpath( datadir, welldir ))
    files = readdir(welldir)
    _, wellname = splitdir(welldir)
    regex = Regex(".*($(escape_regex_str(wellname)).*?)\\d+.mat", "i")
    matches = [match( regex, f ) for f in files]
    prefixes = Set( m[1] for m in matches if m !== nothing )
    length(prefixes) > 1 && error("More than one prefix in dir $welldir. Prefixes:\n$(join(prefixes,"\n"))")
    isempty(prefixes) && error("No mat files found in $welldir matching $regex.")
    first(prefixes)
end

# `datadir` argument kept for backwards compatibility
_filepath( welldir, pr, type, idx ) = "$welldir/$type$pr$(string(idx;pad=4)).mat"
function filepath_f( welldir, datadir = nothing )
    isnothing(datadir) || (welldir = joinpath( datadir, welldir ))
    pr = prefix(welldir)
    (type, idx) -> _filepath( welldir, pr, type, idx )
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

function video_prefix(welldir)
    files = readdir(welldir)
    path, dir = splitdir(welldir)
    dir = chop(dir,tail=2)
    matches = [match( Regex("shape($(escape_regex_str(dir)).*?)\\d+.mp4", "i"), f ) for f in files]
    prefixes = Set( m[1] for m in matches if m !== nothing )
    length(prefixes) > 1 && error("More than one prefix in dir $(dir)")
    isempty(prefixes) ? 
        nothing :
        first(prefixes)
end

# `datadir` argument kept for backwards compatibility
_videopath( welldir, pr, idx ) = "$welldir/shape$pr$(string(idx;pad=4)).mp4"
function videopath_f( welldir, datadir = nothing )
    isnothing(datadir) || (welldir = joinpath( datadir, welldir ))
    pr = video_prefix(welldir)
    pr === nothing && throw("No shape videos in $welldir")
    idx -> _videopath( welldir, pr, idx )
end

read_video(path_f, idx) = VideoIO.load(path_f(idx))

nfiles(path_f) = first( i for i in Iterators.countfrom(0)
                             if !isfile(path_f("coord",i)) )
#                             if !isfile(path_f("corrd",i)) )

using DataFrames
using FileIO

_coords_jld2_path(wellpath) = joinpath(wellpath, "coords_and_size.jld2")
_coords_jld2_path(ex, well, root) = _coords_jld2_path(joinpath(root, ex, well))
function _import_coords_jld2(ex, well, jld2path; err_on_wrong_well)
    @info "Loading coordinates from $jld2path"
    d = load(jld2path)
    if err_on_wrong_well
        @assert ex == d["ex"]
        @assert well == d["well"]
    end
    _ex, _well = get(d, "ex", nothing), get(d, "well", nothing)
    _ex == ex || @warn "Wrong experiment name: expected $(repr(ex)), found $(repr(_ex))."
    _well == well || @warn "Wrong well name: expected $(repr(well)), found $(repr(_well))."
    d["traj"]
end

function load_coords_and_size(ex, well, root; err_on_wrong_well=true)
    jld2path = _coords_jld2_path(ex, well, root)
    if isfile(jld2path)
        return _import_coords_jld2(ex, well, jld2path; err_on_wrong_well)
    else
        @info "No coords file at $jld2path. Loading from mat files."
        df, _ = import_coords(joinpath(ex,well), root; with_size=true)
        return df
    end
end

# Old API, importing directly from `mat` files:
function import_coords( welldir, datadir = nothing; with_size=false )
    path_f = filepath_f(welldir, datadir)
    x = Union{Float64,Missing}[]
    y = Union{Float64,Missing}[]
    with_size && (sz = Union{Float64,Missing}[])
    fileno = Int[]
    batch_boundaries = Int[]
    n = 0
    @progress "Importing .mat files from $welldir" for i in 0:nfiles(path_f)-1
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

boundaries_from_traj(traj) = [searchsortedfirst(traj.fileno, i) - 1 for i in 1:traj.fileno[end]+1]

const frames_per_s = 3

# TODO pixel size per experiment? (pixels were 10μm in older experiments)
# const pixel_size_μm = 10
