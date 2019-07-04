using TOML

const stage_names = ["L1", "L2", "L3", "L4", "A"]
const stages_filepath = normpath("$(@__DIR__)/../stages.toml")
loadstages() = TOML.parsefile(stages_filepath)
savestages(stages) = open( io->TOML.print(io,stages), stages_filepath, "w" )
#appendstages(stages...) = open( io->TOML.print(io,Dict(stages)), stages_filepath, "a" )

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
