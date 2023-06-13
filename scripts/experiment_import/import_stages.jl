"""
Import stages data saved from MATLAB script.
It's better to run `save_to_single_file.jl` first, to make this script faster.
"""

## Find coordinate and analysis output paths

using Elegans
using Elegans: loadstages, getstages!, import_coords
using ScanDir
using Plots
using StatsBase
using DataStructures

function expname(dirname)
    #m = match(r"(RA(\d+)_\d+)", dirname)
    #m = match(r"(NG(\d+))", dirname)
    m = match(r"(SI(\d+))", dirname)
    m === nothing ? nothing : (index = parse(Int, m.captures[2]), name = m.captures[1])
end
function exp_dict(dirs; expname)
    exps = DefaultDict(()->[])
    for d in dirs
        e = expname(d)
        e === nothing && continue
        push!(exps[e.index], (name = e.name, path = d))
    end
    exps
end

# coord_data_root = "U:\\experiments\\reemy"
# analyzed_root = "U:\\experiments\\reemy\\analyzed\\N2_cat2_tph1_050921"
#coord_data_root = "H:/"
# coord_data_root = "U:/experiments/reemy/new"
# analyzed_root = "H:/separate to dev stages/"

# coord_data_root = "U:/experiments/manal/Coords"
# analyzed_root = "U:/experiments/manal/Analyzed"

# coord_data_root = "G:/experiments/sharonin/coords"
# analyzed_root = "G:/experiments/sharonin/roaming"

# coord_data_root = "I:/CB1611, mec-4 experiments/coords"
# analyzed_root = "I:/CB1611, mec-4 experiments/analyzed"


coord_data_root = "G:/experiments/Nabeel/XY Results all"
#analyzed_root = "G:/experiments/Nabeel/separated"
# Nabeel's experiments have a subset of worms re-analyzed in this directory:
#analyzed_root = "G:/experiments/Nabeel/RNA-seq reanalysis for 2 singleand 11 pairs of worms, 21-05-2023"
# Merged with re-analyzed:
analyzed_root = "G:/experiments/Nabeel/separated-merged"

# coord_data_root = "E:\\eshkar"
# analyzed_root = "E:\\eshkar\\analyzed"
# exp_name_re = r"(?:( |_))EN(\d+)"

coord_dirs = readdir(coord_data_root)
#coord_dirs = "NG034", "NG035"
analyzed_dirs = [e.root for e in ScanDir.walkdir(analyzed_root; prune=isfile)]
expname.(coord_dirs)
##
exp_coords = exp_dict(coord_dirs; expname)
exp_analyzed = exp_dict(analyzed_dirs; expname)

all_exps = sort!(collect(union(keys(exp_coords), keys(exp_analyzed))))

## Excel-importable tab-separated table to clipboard

sprint() do io; for ex in all_exps
    c, a = exp_coords[ex], exp_analyzed[ex]
    print(io, "$ex\t$(length(c))\t$(length(a))")
    names = union((x.name for x in c), (x.name for x in a))
    isempty(names) || print(io, "\t", only(names))
    isempty(a) || print(io, "\t", join((x.path for x in a), "\t"))
    println(io)
end; end |> clipboard

## Extract stage boundaries

using ProgressLogging
#using Glob
#bounds_dict = Dict()

to_skip = []

#stagedict = Dict() # loadstages()
stagedict = loadstages()

@progress "exps" for i in setdiff(all_exps, to_skip)
    c, a = exp_coords[i], exp_analyzed[i]
    if isempty(c)
        println("$i: no coord data. skipping.")
        continue
    end
    if isempty(a)
        println("$i: no analyzed data. skipping.")
        continue
    end

    ex_name, ex_path = only(c)
    println("$i: $ex_name")
    c_path = joinpath(coord_data_root, only(c).path)
    #wellpaths = union(glob("CAM*", c_path), glob("*/CAM*", c_path))
    wellnames = [d.name for d in scandir(c_path) if isdir(d) && startswith(d.name, r"CAM"i)] # filter(startswith(r"CAM"i), readdir(c_path))
    isempty(wellnames) && @warn "No wells in $c_path"
    matpaths = [x.path for x in a]
    #@show matpaths
    ex_bounds = Dict()
    # @progress "wells" for wellpath in wellpaths
    #     _, wellname = splitdir(wellpath)
    @progress "wells" for wellname in wellnames
        # wellstages = get( get!( stagedict, ex_name, Dict() ), wellname, nothing )
        # if wellstages !== nothing
        #     @info "Stages found for $ex_name/$wellname. Skipping."
        #     continue
        # end
        well = Well(coord_data_root, ex_path, wellname)
        println("    $wellname")
        # There may be several matpaths. Only one should contain info for this cam.
        relevant_matpaths = filter(matpath->any(fname->startswith(fname, Regex("coord\\Q$wellname\\E","i")),
                                                readdir(matpath)), matpaths)
        if isempty(relevant_matpaths)
            @warn "No analyzed output for $ex_name/$wellname" matpaths
            continue
        elseif length(relevant_matpaths) > 1
            @warn "Multiple mat paths for $ex_name/$wellname: $relevant_matpaths"
        end
        matpath_for_well = first(relevant_matpaths)

        #traj, file_boundaries = import_coords(relpath(wellpath, coord_data_root), coord_data_root)
        traj = load_coords_and_size(well)
        file_boundaries = Elegans.boundaries_from_traj(traj)
        
        stage_boundary_files = Elegans.find_boundary_files(wellname, matpath_for_well)
        stages_mat_format = if isempty(stage_boundary_files)
            @info """No boundary file found for $ex_name/$wellname at $matpath_for_well. 
                     Retrieving from separated coordinates."""
            Elegans.OldMATFormat()
        else
            Elegans.NewMATFormat()
        end
        Elegans.getstages!( ex_name, wellname, matpath_for_well, file_boundaries, stagedict; 
                    mat_fmt = stages_mat_format )
    end
end

## Save to TOML file

Elegans.savestages(stagedict)

##
Elegans.savestages(stagedict, "C:\\Users\\sternlab\\tmp-stages.toml")