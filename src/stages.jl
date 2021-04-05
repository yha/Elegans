using TOML
using AtomicFileWrite
using DataStructures

const stage_names = ["L1", "L2", "L3", "L4", "A"]
const stages_filepath = normpath("$(@__DIR__)/../stages.toml")

sortedstages(stages) = SortedDict(k=>SortedDict(v) for (k,v) in pairs(stages))
loadstages(stagefile=stages_filepath) = TOML.parsefile(stagefile)
savestages(stages, stagefile=stages_filepath) = atomic_write( stagefile ) do io
    TOML.print(io, sortedstages(stages))
end
#appendstages(stages...) = open( io->TOML.print(io,Dict(stages)), stages_filepath, "a" )

function stage_frames( ex, cam, stage; stagedict=loadstages() )
    stage_boundaries = stagedict[ex][cam]
    return stage_boundaries[stage]+1:stage_boundaries[stage+1]
end

nstages( ex, cam; stagedict=loadstages() ) = length(stagedict[ex][cam]) - 1

## Loading of stages saved by MATLAB script (separate_to_developmental_stages.m)

using MAT

struct OldMATFormat end
struct NewMATFormat end

# regex-string multiplication seems to be the only official way to escape a regex currently.
# see https://github.com/JuliaLang/julia/issues/6124
escape_regex(str, flags="") = Regex("", flags) * str

function boundary_file( camname, stages_path )
    boundary_files = filter(fname -> endswith(fname, "stages.mat") && startswith(fname, escape_regex("coord$camname","i")),
                            readdir(stages_path))
    isempty(boundary_files) && error("No stage boundaries file found for cam $camname")
    length(boundary_files) > 1 && error("Multiple ($(length(boundary_files))) stage boundaries file found for cam $camname")
    first(boundary_files)
end

function get_stage_boundaries_new_fmt( camname, stages_path, file_boundaries )
    boundspath = joinpath(stages_path, boundary_file(camname, stages_path))
    stage_boundaries, file_start = matopen(boundspath) do file
        read(file, "boundaries", "file_start")
    end
    file_boundaries[Int(file_start)-1] .+ Int.(stage_boundaries)[1,:]
end

function get_stage_boundaries_old_fmt( camname, mat_path, file_boundaries )
    cam_coord_files = filter(fname->startswith(fname, escape_regex("coord$camname","i")),
                             readdir(mat_path))
    nstages = length(stage_names)
    stage_starts = Vector{Int}(undef, nstages)
    stage_ends = Vector{Int}(undef, nstages)
    for (i, stage_name) in enumerate(stage_names)
        stage_file = only(filter(fname->endswith(fname, "$stage_name.mat"), cam_coord_files))
        mat_data = matopen(joinpath(mat_path,stage_file)) do file
            read(file, "coord_$stage_name")
        end
        start_video, start_frame = Int(mat_data[1,1]), Int(mat_data[1,2])
        # Old scripts sometimes produce extra rows filled with zeros when the
        # adult stage is "too short" (<16hr).
        # Zero in the video number is always a mistake, so find the last
        # non-zero there to locate the actual last frame.
        last_frame_index = findlast( i -> mat_data[i,1] != 0, 1:size(mat_data,1) )
        end_video, end_frame = Int(mat_data[last_frame_index,1]), Int(mat_data[last_frame_index,2])
        stage_starts[i] = file_boundaries[start_video-1] + start_frame
        stage_ends[i] = file_boundaries[end_video-1] + end_frame
    end
    for i=1:nstages-1
        # old format had an overlap of one frame between stages
        @assert stage_ends[i] == stage_starts[i+1]
    end
    stage_boundaries = [stage_starts[1]; stage_ends]
end


get_stage_boundaries(::NewMATFormat, args...) = get_stage_boundaries_new_fmt(args...)
get_stage_boundaries(::OldMATFormat, args...) = get_stage_boundaries_old_fmt(args...)

function getstages!( ex, cam, stages_path, file_boundaries, allstages=loadstages(); mat_fmt=NewMATFormat())
    camstages = get( get!( allstages, ex, Dict() ), cam, nothing )
    if camstages === nothing
        @info "retrieving and storing stages for $ex/$cam from $stages_path"
        camstages = get_stage_boundaries( mat_fmt, cam, stages_path, file_boundaries )
        allstages[ex][cam] = camstages
        #Elegans.savestages(allstages)
    end
    camstages
end


## TODO custom Stages type and AbstractTrees implementation
#using AbstractTrees
#AbstractTrees.children(d::Dict) = collect(d)
#AbstractTrees.children(a::AbstractArray) = ()
