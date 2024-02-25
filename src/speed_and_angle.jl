include("funcs.jl")
using Statistics

function calc_stats!(traj, dframes)
    dt = dframes / frames_per_s
    d_fwd(x) = (res = similar(x);
                res[end-dframes+1:end] .= missing;
                res[begin:end-dframes] .= @views (x[begin+dframes:end] .- x[begin:end-dframes]);
                res)
    d_back(x) = (res = similar(x);
                res[1:dframes] .= missing;
                res[dframes+1:end] .= @views (x[1+dframes:end] .- x[1:end-dframes]);
                res)

    # velocities in pixels/s
    traj.vx = d_fwd(traj.x) / dt
    traj.vy = d_fwd(traj.y) / dt
    traj.speed = [hypot(x,y) for (x,y) in zip(traj.vx, traj.vy)]
    traj.angle = [atan(y,x) for (x,y) in zip(traj.vx, traj.vy)]
    traj.dangle = fix_angle.( d_back(traj.angle) )

    traj
end

function import_and_calc(cam, dframes = 3, datadir = datadir)
    traj, _ = import_coords(cam, datadir)
    calc_stats!(traj, dframes)
end
