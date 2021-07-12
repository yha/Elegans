using Revise
using Elegans
using Images, ImageFiltering
using ProgressLogging
using ScanDir

contour_method = Thresholding(1.0, 0.34)
contours_path = "U:/cached-data/contours"

##

function retry_contours(contours)
    errors_i = sort!(findall( c -> !(c isa AbstractVector), contours.cache ))
    @info "Retrying $(length(errors_i)) error frames..."
    n_errs = 0
    @progress "Retrying contours" for i in errors_i
        delete!(contours.cache, i)
        c = contours(i)
        @assert contours.cache[i] == c
        c isa AbstractVector || (n_errs += 1)
    end
    @info " ... $n_errs errors."
    length(errors_i), n_errs
end


root = "U:/experiments/reemy/"

root = "U:/experiments/reemy"
exs = filter(s->contains(s, r"^RA\d{5}_\d{6}$"), readdir(root))
for ex in exs
    @assert isdir(joinpath(root,ex))
end

exwells  = [(ex, well) for ex in exs for well in readdir(joinpath(root,ex))
                      if isdir(joinpath(root,ex,well)) && startswith(well,r"cam"i)]

##
@progress wells = [Well(root, ex, well) for (ex,well) in exwells]

##
@progress "wells" for well in wells
    println(joinpath(well.experiment, well.well))
    contours, contours_file, vcache = init_contours(well, contour_method, contours_path)
    if !isfile(contours_file)
        @info "No contour file $contours_file. Skipping"
        continue
    end
    errs_before, errs_after = retry_contours(contours)
    errs_before == 0 && continue
    @info "Saving contours for $(well.well) to $contours_file"
    @time save_contours(contours, contours_file)
end

