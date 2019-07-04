include("readdata.jl")
include("funcs.jl")
include("stages.jl")
include("speed_and_angle.jl")

using RollingFunctions

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

function roam_for_stage(cam, traj, stage)
    slope = [5, 2.5, 2.3, 2, 1.5][stage]
    idxs = [round.(Int, stages[cam] .* frames_per_hr); nrow(traj)]
    rows = idxs[2stage-1]:idxs[2stage]
    roam(traj, rows, slope)
end

##

using DataStructures
stages = OrderedDict(sort(pairs(stages)))
n_bins = 75
cams = collect(keys(stages))
fracs = Matrix{Vector{Union{Missing,Float64}}}(undef, length(cams), 5)
for (i,cam) in enumerate(cams)
    @show cam
    traj = import_and_calc(cam, 3)
    for stage = 1:5
        r = roam_for_stage(cam, traj, stage)
        fracs[i,stage] = mean.(Iterators.partition(r,ceil(Int,length(r)/n_bins)))
    end
end

##

using Plots

plot_cams = [1,2,3]
plot_fracs = fracs[plot_cams,:]
plts = [plot(fracs[i,stage], title="$(cams[i]), $(stage_names[stage])", ylims=(0,1))
        for stage in 1:5, i in plot_cams]
plot( plts..., layout=size(plts'), legend=false, titlefontsize=8)

##
frc = reduce(hcat, reduce(vcat,fracs[i,:]) for i=1:length(cams))
heatmap(m2n.(frc'))
vline!((1:4) .* 75, label="", c="white", ls=:dash)

##
hm_roaming(stage) = heatmap( m2n.(reduce(hcat,fracs[:,stage])'),
                             yticks = 1:size(fracs,1),
                             yformatter = x->cams[round(Int,x)],
                             clims = (0,1),
                             title = stage_names[stage] )
using Interact
@manipulate for stage = 1:5
    hm_roaming(stage)
end

##

# function fractions(traj, rows, slope)
#     rm = roam(traj, rows, slope)
#     (roaming = mean(isequal(true),rm),
#      dwelling = mean(isequal(false),rm),
#      missing = mean(ismissing,rm))
# end
#
# function stage_fractions(cam)
#     traj, _ = import_coords("CAM207A2")
#     calc_stats!(traj, 3)
#     n = nrow(traj)
#
#     [fractions(traj, stage_rs(cam,stage,n)...) for stage=1:5]
#         #for (i,slope) in zip(Iterators.partition(idxs,2), slopes)]
# end
#
# ##
# using StatsPlots
# fracs = stage_fractions("CAM207A2")
# plot(stage_names, [f.roaming for f in fracs], label="roaming")
# plot!(stage_names, [f.missing for f in fracs], label="missing")
#
# ##
#
# roam(traj, 10_000:nrow(traj),
