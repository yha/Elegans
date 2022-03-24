include("funcs.jl")
using Statistics

#using DSP

function calc_stats!(traj, dframes)
    dt = dframes / frames_per_s
    #smooth(x) = filtfilt(digitalfilter(Lowpass(1/dframes),FIRWindow(hanning(dframes+2))), x)
    #smooth(x) = filtfilt([0.25,0.5,0.25], x)
    #smooth(x) = filtfilt(fill(1/dframes,dframes), x)
    #smooth(x) = filt(fill(1/dframes,dframes), x)
    #d_pad(x) = [missing; diff(smooth(x)) .* dframes]
    #d(x) = @views (x[1+dframes:end] .- x[1:end-dframes])
    #d_pad(x) = [fill(missing,dframes); d(x)]
    d_fwd(x) = (res = similar(x);
                # res[1:dframes] .= missing;
                # res[dframes+1:end] .= @views (x[1+dframes:end] .- x[1:end-dframes]);
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
    # traj.ax = d_pad(traj.vx) / dt
    # traj.ay = d_pad(traj.vy) / dt
    traj.speed = [hypot(x,y) for (x,y) in zip(traj.vx, traj.vy)]
    # traj.acc = [hypot(x,y) for (x,y) in zip(traj.ax, traj.ay)]
    traj.angle = [atan(y,x) for (x,y) in zip(traj.vx, traj.vy)]
    # traj.dangle = fix_angle.( [missing;diff(traj.angle)] )
    traj.dangle = fix_angle.( d_back(traj.angle) )

    # traj.acc_left = [(-ax*vy + ay*vx)/s for (ax,ay,vx,vy,s) in
    #                                 zip(traj.ax, traj.ay, traj.vx, traj.vy, traj.speed)]
    traj
end

function import_and_calc(cam, dframes = 3, datadir = datadir)
    traj, _ = import_coords(cam, datadir)
    calc_stats!(traj, dframes)
end
