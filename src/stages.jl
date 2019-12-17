using TOML

const stage_names = ["L1", "L2", "L3", "L4", "A"]
const stages_filepath = normpath("$(@__DIR__)/../stages.toml")
loadstages(stagefile=stages_filepath) = TOML.parsefile(stagefile)
savestages(stages, stagefile=stages_filepath) = open( io->TOML.print(io,stages), stagefile, "w" )
#appendstages(stages...) = open( io->TOML.print(io,Dict(stages)), stages_filepath, "a" )


##

## TODO custom Stages type and AbstractTrees implementation
#using AbstractTrees
#AbstractTrees.children(d::Dict) = collect(d)
#AbstractTrees.children(a::AbstractArray) = ()
