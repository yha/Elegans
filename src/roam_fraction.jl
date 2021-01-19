using RollingFunctions

const binsize_speed = 6
const binsize_da = deg2rad(3.6)
const roaming_slope_per_stage = [5, 2.5, 2.3, 2, 1.5]

# TODO fix to smooth in centered windows (no lag), matching MATLAB
_isroaming(speed,dangle,slope) = slope*(speed/binsize_speed + 1) > dangle/binsize_da + 1
function isroaming(traj, rows, slope)
    smooth_window = 10*frames_per_s + 1
    smspeed = runmean(traj.speed[rows], smooth_window)
    smdangle = runmean(abs.(traj.dangle[rows]), smooth_window)
    _isroaming.(smspeed, smdangle, slope)
end

function roam_for_stage(ex, cam, traj, stage; stagedict=loadstages())
    slope = roaming_slope_per_stage[stage]
    idxs = stagedict[ex][cam]
    rows = idxs[stage]:idxs[stage+1]
    isroaming(traj, rows, slope)
end
