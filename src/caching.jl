summarize_exception( e, trace ) = sprint(e, trace) do io, e, trace
    showerror(io, e)
    println(io)
    show(io, MIME("text/plain"), stacktrace(trace))
end

# exc_f defaults to summarize_exception because
#  - BoundError-s on images keep a reference to the image,
#    which makes them very large
#  - Converting such BoundError-s directly to string is extremely
#    slow (~1min!), so not a good option either
function trying_cache(f, ::Type{IN}, ::Type{OUT};
                      cache = Dict{IN,Union{OUT,String}}(),
                      exc_f = summarize_exception) where {IN, OUT}
    x -> get!(cache, x) do
        try
            f(x)
        catch e
            exc_f( e, catch_backtrace() )
        end
    end
end

## Contour caches

using JLD2
using FileIO
using Juno: @progress


# # conversion methods to allow opening older JLD2 files:
# #  - Closed2DCurve.vertices changed from Array to CircularArray
# #  - Exceptions in contour computation now stored as strings, rather than as their types
# #  - Moved from GeomertyTypes.Point to GeometryBasics.Point
# Base.convert(::Type{CircularArray{T,N}}, a::AbstractArray{T,N}) where {T,N} = CircularArray(a)
# Base.convert(::Type{String}, e::Exception) = sprint(showerror, e)
# Base.convert(::Type{Union{Elegans.ContourVec,String}}, e::Exception) = sprint(showerror, e)
# Base.convert(::Type{Union{Elegans.ContourVec,String}}, E::Type{<:Exception}) = string(E)
# Base.convert(::Type{<:CircularVector{<:GeometryBasics.Point{N,T}}},
#              a::CircularVector{GeometryTypes.Point{N,T}}) where {N,T} = CircularVector(GeometryBasics.Point.(a))


#contours_path = "\\\\132.68.111.44\\LabData\\yharel\\contours"
const default_contours_path = joinpath(datadir,"contours")
const ContourVec = Vector{Closed2DCurve{Float64}}

function frame_contours( vcache, i, method )
    fr = centered(get_frame(vcache,i))
    c = raw_worm_contours(Gray.(fr), method)
    #c = Elegans.raw_worm_contours(imfilter(Gray.(fr),Kernel.gaussian(g)), th)
end

contour_cache( videocache, method) = trying_cache(
    i -> frame_contours( videocache, i, method ),
    Int, ContourVec;
    cache = Dict{ Int, Union{ContourVec, String} }(),
    exc_f = summarize_exception )

contours_methodname(m::Thresholding) = "$(m.level)-$(m.σ)"
contours_methodname(m::ContouringMethod) = string(m)

function contours_filename(ex, root, method, contours_path = default_contours_path )
    rel_name = replace( relpath(ex,root),  r"[\\/]" => "-" )
    m = contours_methodname(method)
    contours_file = joinpath(contours_path,"contours-$m-$rel_name.jld2")
end

function load_contours( contours_file, vcache )
    @info "Loading cached contours from $contours_file ..."
    @time d = load(contours_file)
    stored_contours = get(d, "contours", nothing)
    if stored_contours === nothing
        # try older format, which has uncentered contours
        @info "... centering contours ..."
        raw_contours = d["contours.cache"]
        m, n = size(get_frame(vcache, 1))
        midframe = Point2(m÷2+1, n÷2+1)
        # center contours
        @time stored_contours = Dict{Int,valtype(raw_contours)}(
            Base.Generator(raw_contours) do (i,c)
                if c isa AbstractVector # successful frame
                    i => [Closed2DCurve(c.-midframe) for c in raw_contours[i]]
                else
                    i => c
                end
            end)
    end
    n, n_err = length(stored_contours), count((!isa).(values(stored_contours),Vector))
    @info "... $(n-n_err) contours loaded ($n_err errors)"
    stored_contours
end

function init_contours( campath, root, method, contours_path = default_contours_path )
    relcam = replace(relpath(campath,root), "\\"=>"/")
    @info "Initializing video cache ($campath)..."
    vcache = VideoCache(relcam,root)
    contours = contour_cache(vcache,method)

    contours_file = contours_filename( campath, root, method, contours_path )

    if isfile(contours_file)
        stored_contours = load_contours( contours_file, vcache )
        merge!(contours.cache, stored_contours)
    end
    contours, contours_file, vcache
end

save_contours(contours, contours_file) = @save contours_file contours.cache

function compute_all_contours( ex, root, th, g, contours_path = default_contours_path )
    contours, contours_file, vcache = init_contours(ex, root, th, g, contours_path)
    @progress "computing contours" cont = [contours(i) for i = 1:nframes(vcache)]
    contours, contours_file, vcache
end
