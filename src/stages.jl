using TOML

const stage_names = ["L1", "L2", "L3", "L4", "A"]
const stages_filepath = normpath("$(@__DIR__)/../stages.toml")
loadstages(stagefile=stages_filepath) = TOML.parsefile(stagefile)
savestages(stages, stagefile=stages_filepath) = open( io->TOML.print(io,stages), stagefile, "w" )
#appendstages(stages...) = open( io->TOML.print(io,Dict(stages)), stages_filepath, "a" )


## Loading of stages saved by MATLAB script (separate_to_developmental_stages.m)

using MAT

function boundary_file( camname, stages_path )
    boundary_files = filter(fname -> endswith(fname, "stages.mat") && startswith(fname, Regex("coord$camname","i")),
                            readdir(stages_path))
    isempty(boundary_files) && error("No stage boundaries file found for cam $camname")
    length(boundary_files) > 1 && error("Multiple ($(length(boundary_files))) stage boundaries file found for cam $camname")
    first(boundary_files)
end

function get_stage_boundaries( camname, stages_path, file_boundaries )
    boundspath = joinpath(stages_path, boundary_file(camname, stages_path))
    stage_boundaries, file_start = matopen(boundspath) do file
        read(file, "boundaries", "file_start")
    end
    file_boundaries[Int(file_start)-1] .+ Int.(stage_boundaries)[1,:]
end

function getstages!( ex, cam, stages_path, file_boundaries, allstages=loadstages() )
    camstages = get( get!( allstages, ex, Dict() ), cam, nothing )
    if camstages === nothing
        @info "retrieving and storing stages for $ex/$cam from $stages_path"
        camstages = get_stage_boundaries( cam, stages_path, file_boundaries )
        allstages[ex][cam] = camstages
        Elegans.savestages(allstages)
    end
    camstages
end

## TODO custom Stages type and AbstractTrees implementation
#using AbstractTrees
#AbstractTrees.children(d::Dict) = collect(d)
#AbstractTrees.children(a::AbstractArray) = ()
