using ImageFiltering: centered
using Formatting

_try_plot_frame(vcache,i; kw...) = try
    plot(centered(get_frame(vcache,i)); kw...)
catch
    plot(; kw...)
end

plot_abs_frame( vcache, traj, i; kw...) = try
    fr = centered(get_frame(vcache,i))
    plot(traj.x[i] .+ axes(fr,1), traj.y[i] .+ axes(fr,2), fr; kw...)
catch
    plot(; kw...)
end

frametime_str(i, fps::Integer=frames_per_s) = format("{:d}:{:02d}:{:02d}", i÷3600fps, (i÷60fps)%60, (i÷fps)%60)
function plot_frame_with_time(vcache, i; kw...)
    _try_plot_frame(vcache,i; kw...)
    annotate!([(5,65,text(frametime_str(i),:center,:top,"lightgray"))])
end
