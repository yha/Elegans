import TOML
#yield() # workaround Juno.jl bug #316, for running as a Juno cell

include("stages.jl")
include("readdata.jl")
include("speed_and_angle.jl")

using StatsBase: ecdf
self_ecdf(x) = ecdf(collect(skipmissing(x))).(x)

using RollingFunctions

const frames_per_hr = frames_per_s * 3600

function _guess_stages(trace, chunksize, skip, th = 0.2, w = 50)
    min_stage_len_hr = 2
    min_stage_len_chunks = min_stage_len_hr * frames_per_hr / chunksize

    r = runmean(coalesce.(self_ecdf(trace),1),w)
    j = r .< th
    d = findall(diff(j) .!= 0)

    if length(d) >= 9
        b = [d[2]; mean.(partition(d[end-7:end],2)); length(trace) + w/2]
        stage_lens = diff(b)
        minimum(stage_lens) < min_stage_len_chunks && return nothing
        ((b .- w/2) .* chunksize .+ skip)
    else
        nothing
    end
end
function guess_stages(speed, chunksize=300, skip=10_000, th=0.2, w=50)
    trace = mean.(Iterators.partition(speed[skip:end], chunksize))
    _guess_stages(trace, chunksize, skip, th, w)
end

fmtboundaries(boundaries) = join(boundaries ./ frames_per_hr, ", ")
fmtboundaries(::Nothing) = "(none)"

using Interact, Plots
using DataStructures
using ScanDir

const _exclude_dirs = ("\$RECYCLE.BIN", "System Volume Information")

function find_exps(root; onerror=e->(showerror(stderr,e);println(stderr)))
    alldirs = collect(ScanDir.walkdir(root; onerror=onerror,
                                            prune = d -> d.name ∈ _exclude_dirs))
    [relpath(x.root, root) for x in alldirs
            if !isempty(x.dirs) && startswith(lowercase(x.dirs[1]), "cam")]
end

function mark_stages_gui( datadir = datadir, stagefile=stages_filepath )
    stages = DefaultDict(()->Dict{AbstractString,Any}(), loadstages(stagefile))

    hm_toggle = toggle(true, label="hide marked")
    #exps = filter(x -> !(x in _exclude_dirs) && isdir(joinpath(datadir,x)), dir)
    exps = find_exps(datadir)
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
                                    [4,16,24,32,42,58] .* frames_per_hr :
                                    round.(guess ./ frames_per_hr, digits=1) .* frames_per_hr

            t = skip_hrs : chunksize/frames_per_hr : life_hrs
            sel_t = 0:0.1:round(life_hrs,digits=1)
            boundaries = get(stages[exp],cam,nothing)
            new_boundaries = nothing
            btn = button("save")
            on(btn) do _
                TOML.print( stdout, Dict(exp => Dict(cam => new_boundaries)) )
                stages[exp][cam] = new_boundaries
                savestages(stages, stagefile)
                hm_toggle[] = hm_toggle[] # trigger reload
            end
            vbox( fmtboundaries(boundaries), begin
                @manipulate for boundaries_hr_txt = textbox(
                            value = fmtboundaries(something(boundaries, default_boundaries)))
                    plt = plot( t, trace, title = "$exp: $cam",
                                xlabel="hr", ylabel="px/s", legend=false )
                    hr_strings = filter(!isempty, split(boundaries_hr_txt, r"\s*,\s*"))
                    new_boundaries_hr = tryparse.(Float64, hr_strings)
                    if any(isnothing, new_boundaries_hr)
                        new_boundaries = nothing
                        plot!(bg="dimgray")
                    else
                        new_boundaries_fl = clamp.(new_boundaries_hr .* frames_per_hr, 0, nrow(traj))
                        new_boundaries = round.(Int,new_boundaries_fl)
                        @assert maximum(abs.(new_boundaries .- new_boundaries_fl)) < 0.001
                        #@show new_boundaries new_boundaries_hr
                        vline!(new_boundaries_hr)
                    end
                    plot!(plt, t, zeros(length(t)), lw=4, line_z=missings, c=reverse(cgrad(:RdYlGn)) )
                end; end,
                btn
            )
        catch e
            showerror(stderr, e)
            show(stderr, MIME("text/plain"), stacktrace(catch_backtrace()))
            e
        end
    end))
end

using Blink
mark_stages_window(datadir=datadir, stagefile=stages_filepath) = body!(
        Window(), mark_stages_gui(datadir, stagefile))
