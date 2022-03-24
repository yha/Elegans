using ImageFiltering
using LazyArrays

const binsize_speed = 0.9          # = `(speed_max-speed_min) / bins_num_speed` in the MATLAB script
const binsize_da = deg2rad(3.62)   # = `(AV_max-AV_min) / bins_num_AV` in the MATLAB script
const roaming_slope_per_stage = [5, 2.5, 2.3, 2, 1.5]
const runmean_kernel = (len = 10*frames_per_s + 1; centered(fill(1/len,len)))
const runmean_half_len = 5*frames_per_s

#_isroaming(speed,dangle,slope) = slope*(speed/binsize_speed + 1) > dangle/binsize_da + 1

# TODO to match MATLAB scripts:
# - add option to trim 3 frames from stage start and end
# - how to deal with the `speed < max_speed` filtering (which "warps time", changing the effective 
#   time-bin edges) in the notebook?



# `discretize=true` mimics the discretization into bins in the original MATLAB scripts.
# `discretize=false` skips the binning.
_isroaming(speed, dangle, slope; discretize=false) = _isroaming(speed, dangle, slope, 
                                                                discretize ? round : identity)
_isroaming(speed, dangle, slope, f) = slope*(f(speed/binsize_speed) + 1) > f(dangle/binsize_da) + 1

# specify `max_speed` to filter out frames with smoothed speed 0 or higher than max_speed,
# as done in the MATLAB script
function isroaming(traj, rows, slope; max_speed = nothing, discretize = false)
    # smspeed = imfilter(m2n(traj.speed[rows]), runmean_kernel)
    # smdangle = imfilter(m2n(abs.(traj.dangle[rows])), runmean_kernel)
    smspeed = mapwindow(mean, m2n(traj.speed[rows]), -runmean_half_len:runmean_half_len)
    smdangle = mapwindow(mean, m2n(abs.(traj.dangle[rows])), -runmean_half_len:runmean_half_len)
    if max_speed !== nothing
        speed_ok = 0 .< smspeed .< max_speed
        smspeed = smspeed[speed_ok]
        smdangle = smdangle[speed_ok]
    end
    _isroaming.(smspeed, smdangle, slope; discretize)
end

roam_for_stage(well::Well, traj, stage; stagedict=loadstages(), max_speed = nothing, discretize = false) = 
    roam_for_stage(well.experiment, well.well, traj, stage; stagedict, max_speed, discretize)

function roam_for_stage( ex::AbstractString, well::AbstractString, traj, stage; 
                         max_speed = nothing, discretize = false,
                         stagedict=loadstages() )
    slope = roaming_slope_per_stage[stage]
    idxs = stagedict[ex][well]
    rows = idxs[stage]+1:idxs[stage+1]
    isroaming(traj, rows, slope; max_speed, discretize)
end


## Mimicking quirks of the MATLAB code:

# stage frames trimmed as in speed&AV vectors from MATLAB
trimmed_stage_frames(well, stage; stagedict) = stage_frames(well, stage; stagedict)[4:end-3]

# Padding with frames from previous/next stage: the MATLAB code processes the stages glued together
# after each stage is trimmed (although time is not continuous across stages this way)
pre_pad_frames(well, stage; stagedict)  = stage == 1 ? Int[] : trimmed_stage_frames(well, stage-1; stagedict)[end-runmean_half_len+1:end]
post_pad_frames(well, stage; stagedict) = stage == 5 ? Int[] : trimmed_stage_frames(well, stage+1; stagedict)[begin:begin+runmean_half_len-1]

function roam_for_stage_like_matlab( well, traj, stage; 
                                     max_speed = 45,
                                     stagedict=loadstages() )

    fr = ApplyArray(vcat, pre_pad_frames(well, stage; stagedict),
                                 trimmed_stage_frames(well, stage; stagedict),
                                 post_pad_frames(well, stage; stagedict))

    smspeed  = mapwindow(mean, Elegans.m2n(traj.speed[fr]), -runmean_half_len:runmean_half_len)[ begin+runmean_half_len:end-runmean_half_len ]
    smdangle = mapwindow(mean, Elegans.m2n(abs.(traj.dangle[fr])), -runmean_half_len:runmean_half_len)[ begin+runmean_half_len:end-runmean_half_len ]
    speed_ok = 0 .< smspeed .< max_speed
    
    slope = roaming_slope_per_stage[stage]
    roam = _isroaming.(smspeed[speed_ok], smdangle[speed_ok], slope; discretize=true)
    # special edge case for first and last stage: Matlab code leaves the edges as zeros:
    if stage == 1
        roam[1:runmean_half_len] .= 0
    elseif stage == 5
        roam[end-runmean_half_len+1:end] .= 0
    end
           
    (; roam, speed_ok)
end
