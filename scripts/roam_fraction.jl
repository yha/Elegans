using Elegans
using RollingFunctions
using DataFrames

const frames_per_hr = frames_per_s * 3600

const binsize_speed = 6
const binsize_da = deg2rad(3.6)

stages = loadstages()

isroaming(speed,dangle,slope) = slope*(speed/binsize_speed + 1) > dangle/binsize_da + 1

function roam(traj, rows, slope)
    smooth_window = 10*frames_per_s + 1
    smspeed = runmean(traj.speed[rows], smooth_window)
    smdangle = runmean(abs.(traj.dangle[rows]), smooth_window)
    isroaming.(smspeed, smdangle, slope)
end

function roam_for_stage(ex, cam, traj, stage)
    slope = [5, 2.5, 2.3, 2, 1.5][stage]
    idxs = [round.(Int, stages[ex][cam] .* frames_per_hr); nrow(traj)]
    rows = idxs[2stage-1]:idxs[2stage]
    roam(traj, rows, slope)
end

##

using DataStructures
using Statistics
n_bins = 75
exps = sort!(collect(keys(stages)))
cams = mapreduce( pair -> tuple.(pair[1],keys(sort(pair[2]))),
                  append!, sort(stages), init=[])
                  
fracs = Array{Vector{Union{Float64,Missing}},2}(undef,length(cams),5)
@progress for (i,(ex,cam)) in enumerate(cams)
    path = joinpath(ex,cam)
    println(path)
    traj = import_and_calc(path, 3)
    for stage = 1:5
        r = roam_for_stage(ex, cam, traj, stage)
        fracs[i,stage] = mean.(Iterators.partition(r,ceil(Int,length(r)/n_bins)))
    end
end

##

using Plots

function short_cam_name(ex,cam)
    ex_num = match(r"Results(.*)", ex).captures[1]
    cam_num = match(r"CAM(.*)", cam).captures[1]
    "$ex_num: $cam_num"
end

stage_names = Elegans.stage_names

plot_cams = [1,2,3]
plot_fracs = fracs[plot_cams,:]
plts = [plot(fracs[cam,stage],
            title="$(short_cam_name(cams[cam]...)) \n$(stage_names[stage])",
            ylims=(0,1))
        for stage in 1:5, cam in plot_cams]
plot( plts..., layout=size(plts'), legend=false, titlefontsize=8)

##
frc = reduce(hcat, reduce(vcat,fracs[i,:]) for i=1:length(cams))
heatmap(coalesce.(frc',NaN))
vline!((1:4) .* 75, label="", c="white", ls=:dash)

##
hm_roaming(stage) = heatmap( coalesce.(reduce(hcat,fracs[:,stage])', NaN),
                             yticks = 1:size(fracs,1),
                             yformatter = x->short_cam_name(cams[round(Int,x)]...),
                             clims = (0,1),
                             title = stage_names[stage], size=(800,600) )
using Interact
@manipulate for stage = 1:5
    hm_roaming(stage)
end
