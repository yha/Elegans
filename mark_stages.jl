import TOML
yield() # workaround Juno.jl bug #316, for running as a Juno cell

include("stages.jl")
include("readdata.jl")
include("speed_and_angle.jl")
include("funcs.jl")

using StatsBase: ecdf
self_ecdf(x) = ecdf(collect(skipmissing(x))).(x)

using RollingFunctions

function _guess_stages(trace, chunksize, skip, th = 0.2, w = 50)
    frames_per_hr = frames_per_s * 3600
    min_stage_len_hr = 2
    min_stage_len_chunks = min_stage_len_hr * frames_per_hr / chunksize

    r = runmean(coalesce.(self_ecdf(trace),1),w)
    j = r .< th
    d = findall(diff(j) .!= 0)

    if length(d) >= 9
        b = [d[2]; d[end-7:end]]
        stage_lens = diff(b[1:2:end])
        minimum(stage_lens) < min_stage_len_chunks && return nothing
        ((b .- w/2) .* chunksize .+ skip) ./ frames_per_hr
    else
        nothing
    end
end
function guess_stages(speed, chunksize=300, skip=10_000, th=0.2, w=50)
    trace = mean.(Iterators.partition(speed[skip:end], chunksize))
    _guess_stages(trace, chunksize, skip, th, w)
end

fmtboundaries(boundaries) = join(join.(Iterators.partition(boundaries,2),"-"),", ")
fmtboundaries(::Nothing) = "(none)"

using Interact, Plots
using DataStructures

stages = DefaultDict(()->Dict{AbstractString,Any}(), loadstages())

function mark_stages_gui( datadir = datadir )
    stages = DefaultDict(()->Dict{AbstractString,Any}(), loadstages())

    hm_toggle = toggle(true, label="hide marked")
    exps = filter(x->isdir(joinpath(datadir,x)), readdir(datadir))
    exp_dd = dropdown(exps)
    all_cams = Observables.@map filter(fname->isdir(joinpath(datadir,&exp_dd,fname)),
                                       readdir(joinpath(datadir,&exp_dd)))
    cams = Observables.@map(
                &hm_toggle
                    ? filter(x->!haskey(stages[&exp_dd],x), &all_cams)
                    : &all_cams)
    cam_dd = dropdown(cams)

    ui = pad(2em, vbox( hbox(exp_dd, cam_dd, hm_toggle), map(cam_dd) do cam
        try
            exp = exp_dd[]
            traj = import_and_calc(cam, 3, joinpath(datadir,exp_dd[]))
            skip = 10_000
            chunksize = 300
            frames_per_hr = frames_per_s * 3600
            skip_hrs = skip / frames_per_hr
            life_hrs = nrow(traj) / frames_per_hr

            trace = mean.(Iterators.partition(traj.speed[skip:end], chunksize))
            missings = mean.(Iterators.partition(ismissing.(traj.x[skip:end]), chunksize))
            guess = _guess_stages(trace,chunksize,skip)
            default_boundaries = guess === nothing ?
                                    [4,15,17,23,24,31,32,41,42] :
                                    round.(guess,digits=1)

            t = skip_hrs : chunksize/frames_per_hr : life_hrs
            sel_t = 0:0.1:round(life_hrs,digits=1)
            boundaries = get(stages,cam,nothing)
            new_boundaries = nothing
            btn = button("save")
            on(btn) do _
                TOML.print( stdout, Dict(exp => Dict(cam => new_boundaries)) )
                stages[exp][cam] = new_boundaries
                savestages(stages)
                hm_toggle[] = hm_toggle[] # trigger reload
            end
            vbox( fmtboundaries(boundaries), begin
                @manipulate for boundaries_txt = textbox(
                            value = fmtboundaries(something(boundaries, default_boundaries)))
                    plt = plot( t, trace,
                                xlabel="hr", ylabel="um/s", legend=false )
                    m = match(Regex("(.+?)-(.+?),\\s*"^4 * "(.+?)-?\$"), boundaries_txt)
                    if m == nothing
                        new_boundaries = nothing
                        plot!(bg="dimgray")
                    else
                        new_boundaries = tryparse.(Float64, m.captures)
                        new_boundaries = clamp.(something.(new_boundaries, 0), 0, life_hrs)
                        vspan!([new_boundaries;life_hrs], alpha=0.2)
                    end
                    plot!(plt, t, zeros(length(t)), lw=4, line_z=missings, c=:RdYlGn_r )
                end; end,
                btn
            )
        catch e
            e
        end
    end))
end

mark_stages_gui()
