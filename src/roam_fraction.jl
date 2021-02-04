using ImageFiltering

const binsize_speed = 6
const binsize_da = deg2rad(3.6)
const roaming_slope_per_stage = [5, 2.5, 2.3, 2, 1.5]
const runmean_kernel = (len = 10*frames_per_s + 1; centered(fill(1/len,len)))

m2n(v) = replace(v, missing=>NaN)
_isroaming(speed,dangle,slope) = slope*(speed/binsize_speed + 1) > dangle/binsize_da + 1
function isroaming(traj, rows, slope)
    smspeed = imfilter(m2n(traj.speed[rows]), runmean_kernel)
    smdangle = imfilter(m2n(abs.(traj.dangle[rows])), runmean_kernel)
    _isroaming.(smspeed, smdangle, slope)
end

function roam_for_stage(ex, cam, traj, stage; stagedict=loadstages())
    slope = roaming_slope_per_stage[stage]
    idxs = stagedict[ex][cam]
    rows = idxs[stage]+1:idxs[stage+1]
    isroaming(traj, rows, slope)
end
