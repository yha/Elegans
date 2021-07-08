using ImageFiltering

const binsize_speed = 6
const binsize_da = deg2rad(3.6)
const roaming_slope_per_stage = [5, 2.5, 2.3, 2, 1.5]
const runmean_kernel = (len = 10*frames_per_s + 1; centered(fill(1/len,len)))

_isroaming(speed,dangle,slope) = slope*(speed/binsize_speed + 1) > dangle/binsize_da + 1
function isroaming(traj, rows, slope)
    smspeed = imfilter(m2n(traj.speed[rows]), runmean_kernel)
    smdangle = imfilter(m2n(abs.(traj.dangle[rows])), runmean_kernel)
    _isroaming.(smspeed, smdangle, slope)
end

roam_for_stage(well::Well, traj, stage; stagedict=loadstages()) = 
    roam_for_stage(well.experiment, well.well, traj, stage; stagedict)

function roam_for_stage(ex::AbstractString, well::AbstractString, traj, stage; stagedict=loadstages())
    slope = roaming_slope_per_stage[stage]
    idxs = stagedict[ex][well]
    rows = idxs[stage]+1:idxs[stage+1]
    isroaming(traj, rows, slope)
end
