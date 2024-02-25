"""
    Go over experiments and move wells:
        <experiment>/<subdir>/CAMxxxx ==> <expriment>/CAMxxxx
    Also removes the intermediate subdirs (any <experiment>/<subdir> where
    at least one CAM was found).
"""

root = "G:/experiments/Nabeel/XY Results all/"

using Glob
function denest_wells(root; print_only=false)
    for ex in readdir(root)
        expath = joinpath(root,ex)
        println("Reading $expath...")
        #toplevel_campaths = glob("CAM*", expath)
        nested_campaths = glob([fn"*", fn"CAM*"i], expath)

        #isempty(nested_campaths) && continue

        #toplevel_camnames = relpath.(toplevel_campaths,expath)
        intermediates = unique(normpath.(joinpath.(nested_campaths, "..")))

        for campath in nested_campaths
            _, camname = splitdir(campath)
            dest = joinpath(expath, camname)
            println("Move $campath ⟹ $dest")
            print_only || mv(campath,dest)
        end
        for noncam in intermediates
            path = joinpath(expath,noncam)
            # remove .DS_store files too to allow removing otherwise empty dirs
            dsstore = joinpath(path, ".DS_Store")
            if isfile(dsstore)
                println("Remove $dsstore")
                print_only || rm(dsstore)
            end
            println("Remove $path")
            print_only || rm(path)
        end
    end
    #commands
end


# to print commands with moving
denest_wells(root, print_only=true)
##
#denest_wells(root)

