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
    remote, n_local, n_remote
) = read_args(local_project_root)

##
using Elegans

contour_methods = Dict( 1 => Thresholding(1.0,0.35), 
        ( i => Thresholding(1.0,0.34) for i=2:5 )... )

stages = sort!(collect(keys(contour_methods)))

stagedict = loadstages(stage_path)

exs = readdir(local_ex_path)
filter!(in(keys(stagedict)), exs)
@info "Found $(length(exs)) experiments (directories in `$local_ex_path` with entry at `$stage_path`)"

# list of well-ids (experiment, well). These are not Well objects because those store paths which would be different remotely
well_ids = [(ex,well) for ex in exs for well in sort!(collect(keys(stagedict[ex])))]
@info "...$(length(well_ids)) wells"

##
using Distributed

@info "Starting workers..." remote_root remote_project_root local_root local_project_root n_local n_remote

rmprocs(setdiff(workers(), [1])...)
local_procs = addprocs(n_local; exeflags = "--project", dir = local_project_root)
remote_procs = n_remote == 0 ? Int[] :
        addprocs([(remote, n_remote)];
                    dir = remote_project_root, 
                    exename = "$remote_root/julia-$VERSION/bin/julia",
                    exeflags = "--project"
                    )


@info "... added $n_local local, $n_remote remote"


@assert all(==(remote_project_root), remotecall_fetch(pwd, i) for i in remote_procs)

##

# @everywhere using Pkg

# println(join(
#     map(procs()) do i
#         s = remotecall_fetch(()->sprint(io->Pkg.status(;io)), i)
#         "Worker $i:\n$s"
#     end, 
#     "\n"))

# ##

@info "Loading packages..."

@everywhere using Elegans
@everywhere using JLD2
@everywhere using Unzip
@everywhere using DataStructures
@everywhere using ObservablePmap
@everywhere using TerminalLoggers
@everywhere using ProgressLogging
@everywhere using Logging, LoggingExtras


##
@everywhere function maybe_store_midpoints( well, contour_methods, s = Elegans.default_midpoints_s; 
                    midpoints_path, contours_path,
                    stagedict = loadstages(), 
                    headtail_method = Elegans.default_headtail_method, 
                    end_assignment_params = Elegans.EndAssignmentParams() )

    #method2stages = Dict( m => [k for (k,v) in contour_methods if v == m] for m in unique(values(contour_methods)) )
    exname, wellname = well.experiment, well.well
    nstages = Elegans.nstages(well)
    id_str = "$exname/$wellname"
    for method in unique(values(contour_methods))
        stages = [stage for (stage, m) in contour_methods if m == method && stage <= nstages]
        path = midpoints_filename( well.experiment, well.well; contour_method = method,
                                        midpoints_path, headtail_method, end_assignment_params )
        @info "Well $id_str" method stages=string(stages) path
        if isfile(path)
            @info "Midpoints already available" well=id_str method path
        else
            @info "Loading coordinates..." well=id_str
            traj = load_coords_and_size(well)

            @info "Initializing contours..." well=id_str method
            contours, _, _ = init_contours( well, method, contours_path; traj, err_on_missing_contour_file=true )
        
            @info "Computing midpoints..." well=id_str method
            
            midpts, iters, conf = unzip(map(stages) do stage
                @info "Stage $stage (of $stages)" well=id_str method
                irange = stage_frames(well, stage; stagedict)
                range_midpoints( traj, contours, irange, s, headtail_method, end_assignment_params  )
            end)
            stage2midpts = OrderedDict(stages .=> midpts)
            stage2iters = OrderedDict(stages .=> iters)
            stage2conf = OrderedDict(stages .=> conf)

            @info "Storing midpoints..." well=id_str method path
            jldsave(path; midpoints = stage2midpts, iters = stage2iters, conf = stage2conf)
            @info "Midpoints successfully stored" well=id_str
        end
    end
end

# ##
# let (ex, wellname) = well_ids[3]
#     root = local_root
#     ex_root = joinpath(root, ex_dir)
#     well = Well(ex_root, ex, wellname)
#     maybe_store_midpoints(well, contour_methods; stagedict,
#                         contours_path  = "$root/$contours_dir", 
#                         midpoints_path = "$root/$midpoints_dir")
# end

##

using CSSUtil: vbox

using WebIO
using Mux


summ, task = ologpmap(well_ids; logger_f=TerminalLogger, on_error=identity) do (ex, wellname)
    with_logger(EarlyFilteredLogger( log -> log.group != :videoread, current_logger() )) do
        println("$ex: $wellname")
        isremote = myid() ∈ remote_procs
        root = isremote ? remote_root : local_root
        ex_root = joinpath(root, ex_dir)
        @assert isdir(ex_root)

        well = Well(ex_root, ex, wellname)
        maybe_store_midpoints(well, contour_methods; stagedict,
                            contours_path  = "$root/$contours_dir", 
                            midpoints_path = "$root/$midpoints_dir")
        GC.gc()
        println("DONE: $ex: $wellname")
    end
end


h = map(s->HTML("<pre>$s</pre>"), summ)

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
