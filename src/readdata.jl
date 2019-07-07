using MAT
#using Juno
#using ProgressBars

#const datadir = "G:\\Results100616exp\\"
const datadir = normpath("$(@__DIR__)/../data")

function prefix( cameradir, datadir = datadir )
    files = readdir(joinpath( datadir, cameradir ))
    path, dir = splitdir(cameradir)
    prefixes = Set( match( Regex(".*($(dir).*?)\\d+.mat"), f )[1]
                    for f in files )
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
const pixel_size_μm = 10
