"""Save coords, file boundary, and size information to a single file per well
"""

using Elegans
using FileIO
using JLD2

#root = "E:/experiments/reemy"
#exs = filter(s->contains(s, r"^RA\d{5}_\d{6}$"), readdir(root))
root = "U:/experiments/manal/Coords"
exs = filter(s->contains(s, r"^NG\d{3}"), readdir(root))

for ex in exs
    @assert isdir(joinpath(root,ex))
end

exwells  = [(ex, well) for ex in exs for well in readdir(joinpath(root,ex))
                      if isdir(joinpath(root,ex,well)) && startswith(well,r"cam"i)]

## Bad wells:

# badwells = [
#     [
#      ("RA00057_080720", "CAM304A1"),
#      ("RA00031_111219", "CAM302A1"),
#      ("RA00031_111219", "CAM302A3"),
#      ("RA00079_260820", "CAM304A1"),
#      #("RA00013_210719", "CAM035A1"),
#     ];
#     filter(t->t[1] == "RA00012_100619", exwells);   # Two runs (experiment stopped and restarted)
#     filter(t->t[1] == "RA00015_180819", exwells);   # Two runs
#     filter(t->t[1] == "RA00036_300120", exwells)    # Two runs
# ]

badwells = []

#
using ProgressLogging

@progress for (i,(ex,well)) in enumerate(exwells)
    if (ex, well) in badwells
        println( "Skipping bad well $ex / $well" )
        continue
    end
    if isempty(readdir(joinpath(root,ex,well)))
        @warn "Skipping $ex / $well (no files)"
        continue
    end
    dest = joinpath(root, ex, well, "coords_and_size.jld2")
    @assert ispath(dirname(dest))
    if isfile(dest)
        println( "Skipping $ex / $well (already exists: `$dest`)" )
        continue
    end
    println("$i/$(length(exwells)): Loading $ex / $well")
    @time traj, _ = import_coords(joinpath(ex,well), root; with_size=true)
    println("Saving $dest")
    @time @save dest ex well traj
end

