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
                      rethrow_condition = e -> e isa InterruptException,
                      cache = Dict{IN,Union{OUT,String}}(),
                      exc_f = summarize_exception) where {IN, OUT}
    x -> get!(cache, x) do
        try
            f(x)
        catch e
            rethrow_condition(e) && rethrow()
            exc_f( e, catch_backtrace() )
        end
    end
end

## Contour caches

using JLD2
using FileIO
using ProgressLogging: @progress


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


#const default_contours_path = joinpath(datadir,"contours")
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

function contours_filename(well_relpath, method, contours_path )
    wellname = replace( well_relpath,  r"[\\/]" => "-" )
    m = contours_methodname(method)
    contours_file = joinpath(contours_path,"contours-$m-$wellname.jld2")
end
contours_filename(well::Well, method, contours_path) = contours_filename( 
        joinpath(well.experiment, well.well), method, contours_path )

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

function init_contours( well, method, contours_path; traj=load_coords_and_size(well), err_on_missing_contour_file=false )
    @info "Initializing video cache ($(well.path))..."
    vcache = VideoCache(well; traj)
    contours = contour_cache(vcache,method)

    contours_file = contours_filename(well, method, contours_path)

    if isfile(contours_file)
        stored_contours = load_contours( contours_file, vcache )
        merge!(contours.cache, stored_contours)
    elseif err_on_missing_contour_file
        error("No contours file $contours_file")
    end
    contours, contours_file, vcache
end

# TODO use `@save` macro when keyword syntax is supported (next version of JLD2: PR #198):
# @save contours_file contours=contours.cache
save_contours(contours, contours_file) = save(contours_file, Dict("contours"=>contours.cache))

function compute_all_contours( ex, root, th, g, contours_path )
    contours, contours_file, vcache = init_contours(ex, root, th, g, contours_path)
    @progress "computing contours" cont = [contours(i) for i = 1:nframes(vcache)]
    contours, contours_file, vcache
end

## Midpoints cache

#const default_midpoints_path = joinpath(datadir,"midpoints")
const default_midpoints_s = 0:0.025:1

as_tuple(x::T) where T = NamedTuple{fieldnames(T)}(tuple((getfield(x,i) for i in 1:fieldcount(T))...))

const default_headtail_method = SpeedHTCM(5,0)

# `s` argument is unused, but kept for backwards compatibility
# instead of using `s` in the filename, the number of midpoints is inferred
# from the loaded data, with midpoints assumed to be evenly spaced.
function midpoints_filename( ex, wellname, s=nothing; midpoints_path,
        contour_method, headtail_method, end_assignment_params )
    wellname = replace( wellname,  r"[\\/]" => "-" )
    cm = contours_methodname(contour_method)
    m = as_tuple(headtail_method)
    p = as_tuple(end_assignment_params)
    joinpath(midpoints_path,"midpoints-$ex-$wellname $cm $m $p.jld2")
end

save_midpoints(midpts, filename) = save(filename, Dict("midpoints"=>midpts.cache))

function load_midpoints( midpoints_file )
    @info "Loading cached midpoints from $midpoints_file ..."
    @time d = load(midpoints_file)
    stored_midpoints = d["midpoints"]
    n, n_err = length(stored_midpoints), count((!isa).(values(stored_midpoints),Midpoints))
    @info "... midpoints for $(n-n_err) stages loaded ($n_err errors)"
    stored_midpoints
end

# `s` argument is unused, but kept for backwards compatibility
function load_midpoints( well::Well, s=nothing; midpoints_path, contour_method, 
                            headtail_method = default_headtail_method, end_assignment_params=EndAssigmentParams() )
    midpoints_file = midpoints_filename( well.experiment, well.well; midpoints_path, contour_method, headtail_method, end_assignment_params )
    load_midpoints( midpoints_file )
end


const MaybePoint2F = Union{Missing, Point2{Float64}}
const Midpoints = OffsetArray{MaybePoint2F,2,Matrix{MaybePoint2F}}

function midpoint_cache( traj, contours, s = default_midpoints_s;
                headtail_method = default_headtail_method, end_assignment_params = EndAssigmentParams())
    cache = Dict{UnitRange,Midpoints}()
    irange -> get!(cache,irange) do
        midpts, _ = range_midpoints( traj, contours, irange, s, headtail_method, end_assignment_params  )
        midpts
    end
end



# TODO have contour_method stored with contours
# TODO remove this? 
#      a midpoint cache with only some midpoints computed in advance is mostly useless;
#      an error when not all midpoints are pre-computed is better, 
#      so `load_midpoints` and `load_well_midpoints` are more useful.
function init_midpoints( well, traj, contours, s=0:0.025:1;
                        contour_method, midpoints_path,
                        headtail_method = default_headtail_method, end_assignment_params=EndAssigmentParams() )
    mids = midpoint_cache(traj, contours, s; headtail_method, end_assignment_params)

    midpoints_file = midpoints_filename( well.experiment, well.well; midpoints_path, contour_method, headtail_method, end_assignment_params )

    if isfile(midpoints_file)
        stored_midpoints = load_midpoints( midpoints_file )
        merge!(mids.cache, stored_midpoints)
    end
    mids, midpoints_file
end

"""
    load_well_midpoints(ex, well, s, contour_methods, iranges = nothing; 
                        midpoints_path, headtail_method, end_assignment_params = EndAssigmentParams())
Load midpoints for different contouring methods (different files) and merge into a single dict.
`contour_methods` is a dict mapping stage to contouring method.
If `iranges` is given, verify that the loaded frame ranges match the given ranges.
"""
function load_well_midpoints(well, contour_methods, iranges = nothing; 
                        midpoints_path, headtail_method = default_headtail_method, end_assignment_params = EndAssigmentParams())
    method2stages = Dict( m => [k for (k,v) in contour_methods if v == m] for m in unique(values(contour_methods)) )
    mids_dicts = Dict( m => load_midpoints( midpoints_filename( well.experiment, well.well;
                                contour_method = m,
                                midpoints_path, headtail_method,
                                end_assignment_params ) )
                        for m in keys(method2stages) )
    if iranges !== nothing
        for (m,d) in mids_dicts
            loaded = sort!(collect(keys(d)))
            expected = sort!([iranges[k] for k in method2stages[m]])
            loaded == expected || error("Expected ranges $expected, Found $loaded (method $m).")
        end
    end
    reduce(merge, values(mids_dicts))
end