using TOML

const stage_names = ["L1", "L2", "L3", "L4", "A"]

stages_file = "stages.toml"
loadstages() = TOML.parsefile(stages_file)
savestages(stages) = open( io->TOML.print(io,stages), stages_file, "w" )
#appendstages(stages...) = open( io->TOML.print(io,Dict(stages)), stages_file, "a" )

function spans( stages, ex, cam )
    frames_per_hr = frames_per_s * 3600
    n = nrow(worms[ex][cam])
    Iterators.partition(round.(Int, [stages[ex][cam] .* frames_per_hr; n]), 2)
end

##

## TODO custom Stages type and AbstractTrees implementation
#using AbstractTrees
#AbstractTrees.children(d::Dict) = collect(d)
#AbstractTrees.children(a::AbstractArray) = ()
