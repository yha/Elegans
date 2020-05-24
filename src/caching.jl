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

#contours_path = "\\\\132.68.111.44\\LabData\\yharel\\contours"
const default_contours_path = joinpath(datadir,"contours")
const ContourVec = Vector{Closed2DCurve{Float64}}

function frame_contours( vcache, i, g, th )
    fr = get_frame(vcache,i)
    c = Elegans.raw_worm_contours(imfilter(Gray.(fr),Kernel.gaussian(g)), th)
end

contour_cache( videocache, threshold, g ) = trying_cache(
    i -> frame_contours( videocache, i, g, threshold ),
    Int, ContourVec;
    cache = Dict{ Int, Union{ContourVec, String} }(),
    exc_f = summarize_exception )

function contours_filename(ex, root, th, g, contours_path = default_contours_path )
    rel_name = replace( relpath(ex,root),  r"[\\/]" => "-" )
    contours_file = joinpath(contours_path,"contours-$th-$g-$rel_name.jld2")
end

function init_contours( ex, root, th, g, contours_path = default_contours_path )
    relex = replace(relpath(ex,root), "\\"=>"/")
    @info "Initializing video cache ($ex)..."
    vcache = VideoCache(relex,root)
    contours = contour_cache(vcache,th,g)

    #contours_file = joinpath(contours_path,"contours-$th-$g-$(replace(relex,"/"=>"-")).jld2")
    contours_file = contours_filename( ex, root, th, g, contours_path )

    if isfile(contours_file)
        @info "Loading cached contours from $contours_file ..."
        stored_contours = load(contours_file, "contours.cache")
        merge!(contours.cache, stored_contours)
        n, n_err = length(stored_contours), count((!isa).(values(stored_contours),Vector))
        @info "... $(n-n_err) contours loaded ($n_err errors)"
    end
    contours, contours_file, vcache
end

save_contours(contours, contours_file) = @save contours_file contours.cache

function compute_all_contours( ex, root, th, g, contours_path = default_contours_path )
    contours, contours_file, vcache = init_contours(ex, root, th, g, contours_path)
    @progress "computing contours" cont = [contours(i) for i = 1:nframes(vcache)]
    contours, contours_file, vcache
end
