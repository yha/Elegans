using TOML

const stage_names = ["L1", "L2", "L3", "L4", "A"]

stages_file = "stages.toml"
loadstages() = TOML.parsefile(stages_file)
savestages(stages) = open( io->TOML.print(io,stages), stages_file, "w" )
appendstages(stages...) = open( io->TOML.print(io,Dict(stages)), stages_file, "a" )
