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

# coord_data_root = "E:\\eshkar"
# analyzed_root = "E:\\eshkar\\analyzed"
coord_data_root = "U:\\experiments\\eshkar"
analyzed_root = "U:\\experiments\\eshkar\\analyzed"

function expname(dirname)
    # # match two uppercase letters followed by a sequence of digits
    # m = match(r"(\p{Lu}{2}(\d+))", dirname)
    #m === nothing ? nothing : (index = parse(Int, m.captures[2]), name = m.captures[1])
    m = match(r"(EN\d+)", dirname)
    m === nothing ? nothing : m.captures[1]
end


function exp_dict(dirs; expname)
    exps = DefaultDict(()->String[])
    for d in dirs
        e = expname(d)
        e === nothing && continue
        #push!(exps[e.index], (name = e.name, path = d))
        push!(exps[e], d)
    end
    exps
end


coord_dirs = readdir(coord_data_root)
analyzed_dirs = [e.root for e in ScanDir.walkdir(analyzed_root; prune=isfile)]
expname.(coord_dirs)
##
exp_coords = exp_dict(coord_dirs; expname)
exp_analyzed = exp_dict(analyzed_dirs; expname)
@assert keys(exp_analyzed) ⊆ keys(exp_coords)

for (k,v) in exp_coords
    length(v) > 1 && error("More than one coords directory for experiment $k: $v")
end

all_exps = sort!(collect(union(keys(exp_coords), keys(exp_analyzed))))

## Excel-importable tab-separated table to clipboard

sprint() do io; for ex in all_exps
    c, a = exp_coords[ex], exp_analyzed[ex]
    print(io, "$ex\t$(length(c))\t$(length(a))")
    isempty(c) || print(io, "\t", only(c))
    print(io, "\t", join(a, "\t"))
    println(io)
end; end |> clipboard

## Extract stage boundaries

using ProgressLogging

to_skip = []

stagedict = loadstages()

@progress "exps" for ex_name in setdiff(all_exps, to_skip)
    coords, analyzed = exp_coords[ex_name], exp_analyzed[ex_name]
    if isempty(coords)
        println("$ex_name: no coord data. skipping.")
        continue
    end
    if isempty(analyzed)
        println("$ex_name: no analyzed data. skipping.")
        continue
    end

    ex_path = only(coords)
    println("$ex_name @ $ex_path")
    c_path = joinpath(coord_data_root, only(coords))
    wellnames = [d.name for d in scandir(c_path) if isdir(d) && startswith(d.name, r"CAM"i)] 
    isempty(wellnames) && @warn "No wells in $c_path"
    #@show analyzed
    ex_bounds = Dict()
    @progress "wells" for wellname in wellnames
        well = Well(coord_data_root, ex_path, wellname)
        println("    $wellname")
        # There may be several analyzed paths. Only one should contain info for this cam.
        relevant_a_paths = filter(analyzed) do path
            any(fname->startswith(fname, r"coord"i * wellname), readdir(path))
        end
        if isempty(relevant_a_paths)
            @warn "No analyzed output for $ex_name/$wellname" analyzed
            continue
        elseif length(relevant_a_paths) > 1
            @warn "Multiple mat paths for $ex_name/$wellname: $relevant_a_paths"
        end
        matpath_for_well = first(relevant_a_paths)

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