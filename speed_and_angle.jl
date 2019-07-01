include("funcs.jl")
using Statistics

using DSP

function calc_stats!(traj, dframes)
    dt = dframes / frames_per_s
    #smooth(x) = filtfilt(digitalfilter(Lowpass(1/dframes),FIRWindow(hanning(dframes+2))), x)
    #smooth(x) = filtfilt([0.25,0.5,0.25], x)
    #smooth(x) = filtfilt(fill(1/dframes,dframes), x)
    #smooth(x) = filt(fill(1/dframes,dframes), x)
    #d_pad(x) = [missing; diff(smooth(x)) .* dframes]
    #d(x) = @views (x[1+dframes:end] .- x[1:end-dframes])
    #d_pad(x) = [fill(missing,dframes); d(x)]
    d_pad(x) = (res = similar(x);
                res[1:dframes] .= missing;
                res[dframes+1:end] .= @views (x[1+dframes:end] .- x[1:end-dframes]);
                res)

    traj.vx = d_pad(traj.x) * pixel_size_μm / dt
    traj.vy = d_pad(traj.y) * pixel_size_μm / dt
    #traj.ax = d_pad(traj.vx) * pixel_size_μm / dt
    #traj.ay = d_pad(traj.vy) * pixel_size_μm / dt
    traj.speed = [hypot(x,y) for (x,y) in zip(traj.vx, traj.vy)]
    #traj.acc = [hypot(x,y) for (x,y) in zip(traj.ax, traj.ay)]
    traj.angle = [atan(y,x) for (x,y) in zip(traj.vx, traj.vy)]
    #traj.dangle = fix_angle.( [missing;diff(traj.angle)] )
    traj.dangle = fix_angle.( d_pad(traj.angle) )

    # traj.acc_left = [(-ax*vy + ay*vx)/s for (ax,ay,vx,vy,s) in
    #                                 zip(traj.ax, traj.ay, traj.vx, traj.vy, traj.speed)]
    traj
end

function import_and_calc(cam, dframes = 3, datadir = datadir)
    traj, _ = import_coords(cam, datadir)
    calc_stats!(traj, dframes)
end
