using Pkg

local_project_root = @__DIR__

Pkg.activate(local_project_root)
@info "Installing packages..."
Pkg.instantiate()

##

include("read_args.jl")
(;
    stage_path, local_root, remote_root, remote_project_root, logdir,
    ex_dir, contours_dir, midpoints_dir, local_ex_path, local_contours_path,
    remote, n_remote, n_local,
    skip_existing
) = read_args(local_project_root)

##

using Elegans
using Distributed


contour_methods = Dict( 1 => Thresholding(1.0,0.34), 
    (i => Thresholding(1.0,0.35) for i in 2:5)...
)

contour_method2stages = Dict(
    (method => sort([stage for (stage,m) in contour_methods if m == method])) for method in unique(values(contour_methods))
)

allstages = loadstages(stage_path)

##

exs = readdir(local_ex_path)
filter!(in(keys(allstages)), exs)
@info "Found $(length(exs)) experiments (directories in `$local_ex_path` with entry at `$stage_path`)"

well_ids = [(ex,well) for ex in exs for well in sort!(collect(keys(allstages[ex])))]
@info "...$(length(well_ids)) wells"

##

using Distributed


@info "Starting workers..." remote_root remote_project_root local_root local_project_root n_local n_remote

rmprocs(setdiff(workers(), [1])...)
local_procs = cd(local_project_root) do
    addprocs(n_local; exeflags = "--project")
end
remote_procs = n_remote == 0 ? Int[] :
                addprocs([(remote, n_remote)];
                      dir = remote_project_root, 
                      exename = "$remote_root/julia-$VERSION/bin/julia",
                      exeflags = "--project"
                      )


@assert all(==(remote_project_root), remotecall_fetch(pwd, i) for i in remote_procs)

# ##
# @everywhere using Pkg

# println(join(
#     map(procs()) do i
#         s = remotecall_fetch(()->sprint(io->Pkg.status(;io)), i)
#         "Worker $i:\n$s"
#     end, 
#     "\n"))


##

@everywhere function contour_counts(well_data, nbins)
    n = length(well_data.traj.x)
    bin_edges = round.(Int, range(0, n; length=nbins+1))
    completed = sort!(collect(keys(well_data.contours.cache)))
    cumcounts = [searchsortedfirst( completed, edge ) for edge in bin_edges]
    counts = diff(cumcounts)
    counts, bin_edges
end


##

@everywhere using ObservablePmap

##

@everywhere using IterTools
@everywhere begin
    block_repr(p, empty="\u2591") = (r=round(Int,8p)) >= 8 ? "█" : r <= 0 ? empty : '\u2580' + r
    _stage_str(len, name) = let s = max(len-1-length(name),0)
        #"\\" * "⋅"^(s÷2) * name * "⋅"^((s+1)÷2) * "/"
        " "^(s÷2) * name * " "^((s+1)÷2)
    end
    function cam_progress_str( well_data, nbins, stage_ends, frames )
        counts, bin_edges = contour_counts(well_data, nbins)
        binlens = diff(bin_edges)
        #rep(p) = (r=round(Int,10p)) >= 10 ? "#" : string(r)
        bin_required = [!isempty(intersect(frames,bin)) for bin in IterTools.partition(bin_edges,2,1)]
        progline = join((block_repr( c/l, r ? " " : "\u2591" )
                        for (c,l,r) in zip(counts,binlens,bin_required)))

        n = length(well_data.traj.x)
        stage_edge_bins = round.(Int, stage_ends .* (nbins/n))
        #@show stage_ends stage_edge_bins
        stageline = join((_stage_str(d,n) for (d,n) in zip(
                        diff([0;stage_edge_bins]), [" "; Elegans.stage_names])), "|") * "|"
        progline * "\n" * stageline * "\n"
    end
end

##


using Elegans
@everywhere using Elegans

##

@everywhere using DataStructures
@everywhere using ProgressLogging
@everywhere using TerminalLoggers

@everywhere function store_contours_remote(well_data, stage_ends, frames; chunklen = 1000)
    set_status(str) = @info str

    computed = collect(keys(well_data.contours.cache))
    todo = setdiff!(OrderedSet(frames), computed)

    set_status("Computing $(length(todo)) out of $(length(frames)) contours.")

    for (i,c) in enumerate(todo)
        well_data.contours(c)
        if mod(i,chunklen) == 0
            str = cam_progress_str(well_data, 100, stage_ends, frames)
            set_status("Contours progress: \n" * str)
        end
    end
    if isempty(todo)
        set_status("No new contours computed.")
    else
        set_status("Saving contours to $(well_data.contours_file)...")
        save_contours(well_data.contours, well_data.contours_file)
        set_status("Done.")
    end
end


##

well_method_combinations = Iterators.product(well_ids, keys(contour_method2stages))

using WebIO, CSSUtil, Mux
@everywhere using Logging: current_logger, with_logger
@everywhere using LoggingExtras: EarlyFilteredLogger

summ, task = ologpmap(well_method_combinations; on_error=identity) do ((ex, wellname), contour_method)
    with_logger(EarlyFilteredLogger( log -> log.group != :videoread, current_logger() )) do

        isremote = myid() in remote_procs
        root = isremote ? remote_root : local_root

        exps_path = joinpath( root, ex_dir )
        contours_path = joinpath( root, contours_dir )
        #@show exps_path contours_path
        @assert isdir(joinpath(exps_path,ex))
        well = Well(exps_path, ex, wellname)
        well_str = "$(well.experiment): $(well.well)"
    
        @info "Loading $well_str..."
        if skip_existing
            contours_file = Elegans.contours_filename( well, contour_method, contours_path )
            if isfile(contours_file)
                @info "Skipping existing contour file: $contours_file"
                return
            end
        end

        traj = load_coords_and_size(well)
        contours, contours_file, vcache = init_contours( well, contour_method, contours_path; traj )
        well_data = (; contours, contours_file, traj)

        @info "... loaded $well_str."
        stage_ends = allstages[well.experiment][well.well]
        stage_ranges = [x+1:y for (x,y) in IterTools.partition(stage_ends,2,1)]
        stages = contour_method2stages[contour_method]
        stage_i = intersect(stages, 1:length(stage_ends)-1)
        frames = reduce(union, stage_ranges[stage_i])
        store_contours_remote( well_data, stage_ends, frames )
    end
end

h = map(x -> style(HTML("<pre>$x</pre>")), summ)

port = rand(8100:8200)
WebIO.webio_serve(page("/", vbox(h)), port)
address = "localhost:$port"
println()
println("Progress available at $address (copied to clipboard)")
clipboard(address)

Threads.@spawn begin
    wait(task)
    @info "Done"
end

nothing

