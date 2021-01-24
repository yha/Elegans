using GLMakie
using ProgressLogging
using ProgressLogging: @progress
using Colors
using GeometryBasics

using Elegans

##
using Missings


# convert Plots.jl-style columns to NaN-separated vectors of points
# (init is a single NaN due to Makie issue #533)
_nan(x) = NaN
_nan(::Type{<:Point{N}}) where N = Point(ntuple(i->NaN,N)...)
nan_separated_cols(a; nan = _nan(nonmissingtype(eltype(a)))) = reduce( vcat,
            Iterators.flatten(zip(eachcol(a), Iterators.repeated([nan])));
            #init=Float64[])
            init = [nan] )
nan_separated_cols(x,y) = reduce( vcat,
            Iterators.flatten(zip(eachcol(Point2.(x,y)), Iterators.repeated([Point2(NaN,NaN)])));
            #init=Point2{Float64}[])
            init = [Point2(NaN,NaN)] )

##
# Synced recipes, to allow atomic update input arguments, preventing double
# update through different listeners.
# `Any`-typed fields to work around Makie bug #804
struct HMData
    x
    y
    z
end
 AbstractPlotting.convert_arguments(P::Type{<:Heatmap}, d::HMData) = convert_arguments(P, d.x, d.y, d.z)

@recipe(fig->Theme(), CLines, cl)
function AbstractPlotting.plot!(plt::CLines)
    cl = plt[:cl]
    lines!( plt, lift(cl->cl.p, cl),
                 color = lift(cl->cl.c, cl),
                 colormap = lift(cl->cl.cmap, cl) )
    plt
end

@recipe(fig->Theme(), SyncImage, data)
function AbstractPlotting.plot!(plt::SyncImage)
    data = plt[:data]
    image!( plt, lift(d->d.x, data),
                 lift(d->d.y, data),
                 lift(d->d.img, data) )
    plt
end

##
using StatsBase
winlen, dwin, s = 60, 20, 0:0.025:1
s_edges = [first(s)-step(s)/2; midpoints(s); last(s)+step(s)/2]

##
contours_path = "U:/cached-data/contours"
midpoints_path = "U:/cached-data/midpoints"


# root = "E:/eshkar"
# ex = "Results_090220_EN00021_2nd_new"
root, ex = "U:/experiments/reemy", "230120_RA00035"
stagedict = loadstages()

cams = sort!(collect(keys(stagedict[ex])))[1:8]

exname = match(r"(^|_)\d{6}_(\w\w\d{5})($|_)", ex).captures[2]
#_vcache = VideoCache(joinpath(ex,cams[1]),root)
#get_frame(_vcache, 700000)

##

# TODO: more pdf estimation methods (ASH, KDE)
histogram_pdf(r,edges) = normalize(fit(Histogram, filter(isfinite,r), edges), mode=:pdf)

function make_pdf_matrix( y, edges, pdf_estimate = histogram_pdf )
    hists = [pdf_estimate(r,edges) for r in eachrow(y)]
    z = reduce(hcat, h.weights for h in hists)'
end

##


using LinearAlgebra
using Unzip

# TODO reorganize to deal with different contour method for L1
function make_cam_data(cam, stages_i;
                    midpoints_path, contour_method, headtail_method, stagedict)
    # @show cam
    iranges = [stage_frames(ex,cam,stage_i;stagedict) for stage_i in stages_i]
    mids_dict = Elegans.load_midpoints( Elegans.midpoints_filename( ex, cam, s;
                                    midpoints_path, contour_method, headtail_method,
                                    end_assignment_params=Elegans.EndAssigmentParams() ) )
    # @show sort!(collect(keys(mids_dict)))
    mids = [mids_dict[irange] for irange in iranges]
    @progress "covs" wncovs = [normed_midpoint_covs( m, s, irange; winlen, dwin )
                                 for (m, irange) in zip(mids, iranges)]
    ncov, windows = unzip( wncovs )
    log10gv = [log10.(det.(n))' for n in ncov]
    # TODO normalize normal speed
    @progress "normal" nspeed = [normal_speeds(m)' for m in mids]

    traj = import_and_calc(joinpath(ex, cam), 3, root)
    vcache = VideoCache(joinpath(ex,cam), root)
    (; traj, iranges, windows, ncov, nspeed, log10gv, mids, vcache )
end

stages_i = 2:5
contour_method = Thresholding(1.0,0.34)
headtail_method = Elegans.SpeedHTCM(5,0)

@progress "cam data" allcams_data = [make_cam_data(cam, stages_i; stagedict, midpoints_path, contour_method, headtail_method)
                for cam in cams]

##

using OffsetArrays: no_offset_view

function latest_change(obs::Observable{T}, eq = (==)) where T
    out = Observable{T}(obs[])
    on(obs) do val
        eq(val, out[]) || (out[] = val)
    end
    out
end

function translated_frameindex(i, source_iranges, dest_iranges)
    stage_i = searchsortedfirst(last.(source_iranges), i)
    source_start, source_end = extrema(source_iranges[stage_i])
    dest_start,   dest_end   = extrema(dest_iranges[stage_i])

    frac = (i - source_start) / (source_end - source_start)
    round(Int, dest_start + frac * (dest_end - dest_start))
end

# Translate the frame range `irange` to equivalent (time-normalized) ranges for
# other cameras
function translated_frameranges(allcams_data, irange, cam_i)
    source_ranges = allcams_data[cam_i].iranges
    dest_ranges = [vd.iranges for vd in allcams_data]
    first_i, last_i = first(irange), last(irange)
    translated_ranges = [k == cam_i ? irange :
                                      UnitRange( translated_frameindex( first_i, source_ranges, d ),
                                                 translated_frameindex( last_i,  source_ranges, d ) )
                         for (k,d) in enumerate(dest_ranges)]
end

# window indices `j` inside given frame range, where the line (`s`,`y[:,j]`)
# passes through the given `rect`
# (meaning `(s[i],y[i,j])` is in `rect` for some `i`)
function filtered_windows(s, windows, y, rect=nothing, framerange=nothing)
    k = Set(eachindex(windows))
    framerange !== nothing && filter!( j -> issubset(windows[j], framerange), k )
    rect !== nothing && filter!( j -> any((pt ∈ rect) for pt in Point2.(s,y[:,j])), k )
    collect(k)
end

# Window indices for reference statistics, for the given rectangle and time range.
# These are indices of windows for each cam in the same stage,
# within the same normalized time as `framerange` in camera `cam_i` time,
# and where the measure `y` lies inside `rect`.
# The given function `f` returns `(windows, y) = f(camdata, stage_i)`
function ref_samples( allcams_data, s, f, cam_i, stage_i, framerange, rect = nothing )
    @assert issubset( framerange, allcams_data[cam_i].iranges[stage_i] )
    ranges = translated_frameranges(allcams_data, framerange, cam_i)
    [filtered_windows(s, f(camdata,stage_i)..., rect, r)
            for (camdata, r) in zip(allcams_data, ranges)]
end

function ref_pdf_estimates( allcams_data, s, f, cam_i, stage_i, framerange, rect = nothing;
                            edges, pdf_estimate = histogram_pdf )
#    @show (f, cam_i, stage_i, framerange, rect, edges)
    @info "Generating reference pdf estimates for $(cams[cam_i]) / $(stages_i[stage_i]) / $framerange..."
    @time samples = ref_samples(allcams_data, s, f, cam_i, stage_i, framerange, rect)
    @progress res = [make_pdf_matrix(f(cd,stage_i)[2][:,k], edges, pdf_estimate) for (cd,k) in zip(allcams_data, samples)]
    @info "Done."
    res
end

# normal speed is estimated from a 3-window around each frame.
function nspeed_windows(cd, stage_i)
    r = cd.iranges[stage_i]
    map( i -> (i-1:i+1) ∩ r, r )
end

log10gv(cd, stage_i) = cd.log10gv[stage_i]
nspeed(cd, stage_i) = no_offset_view(cd.nspeed[stage_i])
log10nspeed(cd, stage_i) = log10.(abs.(no_offset_view(cd.nspeed[stage_i])))

log10gv_visdata(cd, stage_i) = (cd.windows[stage_i], cd.log10gv[stage_i])
nspeed_visdata(cd, stage_i) = (nspeed_windows(cd, stage_i), nspeed(cd, stage_i))
log10nspeed_visdata(cd, stage_i) = (nspeed_windows(cd, stage_i), log10nspeed(cd, stage_i))

##

function make_stage_vis_data(cam_i, stage_i, f, cam_data=allcams_data[cam_i])
    irange = cam_data.iranges[stage_i]

    @info "Generating visualization data for stage $(stages_i[stage_i])..."

    vcache = cam_data.vcache
    traj = cam_data.traj
    sm_speed = imfilter( replace(traj.speed[irange], missing=>NaN), gaussian(100), NA() )
    windows, y = f(cam_data, stage_i)
    @info "Done."

    #miny, maxy = floor(Int,NaNMath.minimum(y)), ceil(Int,NaNMath.maximum(y))
    yf = filter(isfinite, y)
    miny, maxy = floor(Int,quantile(yf,0.005)), ceil(Int,quantile(yf,0.995))
    y_edges = range(miny, maxy; length=41)
    default_rect = FRect2D( 0.0, miny, 1.0, maxy-miny )

    mids = cam_data.mids[stage_i]

    # inputs (cam_i, stage_i, f) included in output to allow synchronized
    # updates through this output, avoiding diamond-shaped dependencies
    (; cam_i, stage_i, f,
        y, mids, traj, sm_speed, windows, irange,
        y_edges, default_rect,
        vcache)
end

using Statistics
using StatsBase
using StatsBase: normalize
using ImageFiltering
using ImageFiltering.KernelFactors: gaussian
using AbstractPlotting.MakieLayout
using Observables
using NaNMath
using GeometryBasics: decompose
using Formatting
using OffsetArrays: no_offset_view

fig = Figure()
# Top menu moved to [1,2] to avoid trigerring hover events on speed plot
# when selecting menu items.
# Can be moved to the top (above speed plot) when when Makie.jl issue #777 is fixed.
top_menus = fig[1,2]

top_menus[1,1] = cam_menu = Menu(fig; options=cams)
top_menus[1,2] = stage_menu = Menu(fig; options=stages_i)
top_menus[1,3] = stats_type_menu = Menu(fig; options=["GV", "normal"], i_selected=1)

cam_menu.i_selected[] = 1
stage_menu.i_selected[] = 1

vis_data_f = lift(i->(log10gv_visdata, log10nspeed_visdata)[i], stats_type_menu.i_selected; typ=Any)

vis_data = lift( make_stage_vis_data,
                 cam_menu.i_selected, stage_menu.i_selected, vis_data_f;
                 typ=Any )

top_menus[1,4] = hm_type_menu = Menu(fig; i_selected=1,
            options = ["x", "x - x̄", "log₂(x/x̄)", "x̄"])
hm_types = (; this=1, Δ=2, logratio=3, mean=4)

hist_time_filter = Toggle(fig; active=true)
hist_sel_filter = Toggle(fig)

fig[1:2,1] = speed_plot = Axis(fig)
hist_area = fig[3,1]

hist_area[1,1] = hist1d_ax = Axis(fig)
hist_area[1,2] = ax = Axis(fig)
hist_area[1,3] = cb_layout = GridLayout()


rowsize!(fig.layout, 2, Auto(true,1))
rowsize!(fig.layout, 3, Auto(true,5))

hist_layout = content(hist_area)
colsize!(hist_layout, 1, Auto(true,1))
colsize!(hist_layout, 2, Auto(true,7))
colsize!(hist_layout, 3, Auto(true,0.25))

linkyaxes!(hist1d_ax, ax)
hideydecorations!(ax)
hidexdecorations!(hist1d_ax)
hidespines!(hist1d_ax)

hist1d_ax.xreversed = true

lines!(speed_plot, lift(d->Point2.(d.irange,d.sm_speed), vis_data))
deregister_interaction!(speed_plot, :rectanglezoom)
function time_segment(irange, from, to)
    i0, i1 = extrema(irange)
    f(t) = i0 + (i1-i0)*t
    round(Int, f(from)):round(Int, f(to))
end


function _fix_range(from, to, irange)
    from, to = round.(Int, minmax(from, to))
    (from:to) ∩ vis_data[].irange
end
frame_to_stagetime(i, irange) = (i-first(irange)) / (last(irange)-first(irange))
stagetime_to_frame(x, irange) = first(irange) + x * (last(irange)-first(irange))
stage_timespan_to_framerange(span, irange) = _fix_range(
                stagetime_to_frame(span[1], irange),
                stagetime_to_frame(span[2], irange),
                irange)

time_sel = select_line(speed_plot.scene)

speed_lims = @lift extrema(filter(isfinite,$vis_data.sm_speed))
#reset_time_sel(sel=time_sel) = (sel[] = Point2.( [first(vis_data[].irange),last(vis_data[].irange)], speed_lims[] ))
reset_time_sel(vd, sel=time_sel) = (s = collect(Point2.( extrema(vd.irange), speed_lims[] ));
                                    #@show s;
                                    sel[] = s)
reset_time_sel(vis_data[])

time_range = lift(time_sel) do line
    @show line
    _fix_range( line[1][1], line[2][1], vis_data[].irange )
end
# Time span as a tuple (`from`,`to`) in normalized stage time
# (0 is stage start, 1 stage end).
# This is kept constant and actual frame range adjusted when switching stages/cams
stage_time_span = lift(time_range) do r
    vd = vis_data[]
    span = (frame_to_stagetime(first(r), vd.irange), frame_to_stagetime(last(r), vd.irange))
    # @show span
    # span
end
on(vis_data) do vd
    # adjust selected frame range according to current normalized stage time
    span = stage_time_span[]
    time_sel[] = collect(Point2.( stagetime_to_frame.(span, Ref(vd.irange)), speed_lims[] ))
end

time_rect = @lift FRect( first($time_range), $speed_lims[1],
                         last($time_range) - first($time_range), $speed_lims[2] - $speed_lims[1] )
poly!(speed_plot, time_rect, color=RGBAf0(0.5,0,1,0.1))
# windows included in selected time range
time_range_i = @lift findall( w -> issubset(w, $time_range), $vis_data.windows )

n_on_text = Textbox(fig; stored_string="100")
n_off_text = Textbox(fig; stored_string="100")
n_on = lift(s->something(tryparse(Int, s), 0), n_on_text.stored_string)
n_off = lift(s->something(tryparse(Int, s), 0), n_off_text.stored_string)
sp_options = fig[2,2] = GridLayout(tellheight=false)
sp_options[1,1] = n_segments_sl = Slider(fig; range=1:100, value=1)
sp_options[2,1] = segment_i_sl = Slider(fig; range=1:1, value=1)
on(n_segments_sl.value) do n
    segment_i_sl.range[] = 1:n
    #set_close_to!(segments_i_sl, segment_i_sl.value[])
    # new_i = min(segment_i_sl.value[], n)
    # segment_i_sl.value[] =  new_i
end
onany(n_segments_sl.value, segment_i_sl.value) do n,i
    frames = time_segment(vis_data[].irange, (i-1)/n, i/n)
    time_sel[] = collect(Point2.( extrema(frames), speed_lims[] ))
end

fig[3,2] = GridLayout(tellheight=false)
rightpane = fig[3,2]

cb_layout[1,1:2] = crange_toggle = Toggle(fig)
cb_layout[2,2] = crange_sl = Slider(fig; range=0:0.1:10,    horizontal=false)

rightpane[1,1:3] = window_choose = GridLayout()
rightpane[2,1:3] = window_viz = GridLayout()
hist_filters_grid = grid!([hist_time_filter Label(fig, "time");
                           hist_sel_filter Label(fig, "rect.")])
on_α, off_α = [Slider(fig; range=0:0.01:1) for _=1:2]
on_α_label, off_α_label = [Label(fig, @lift(format("{:.2f}", $(sl.value))); halign = :right) for sl in (on_α, off_α)]

onoffgrid = grid!([Label(fig, "on")  n_on_text  on_α   on_α_label;
                   Label(fig, "off") n_off_text off_α  off_α_label],
                   width=300)
set_close_to!(on_α, 0.3)
rightpane[3,1] = hist_filters_grid
rightpane[3,2] = onoffgrid
rightpane[3,3] = resample_btn = Button(fig; label="resample")
on(resample_btn.clicks) do _
    active_rect[] = active_rect[]
end

window_choose[1,1] = btn_prev_on = Button(fig; label="⟨")
window_choose[1,2] = btn_prev = Button(fig; label="-")
window_choose[1,3] = window_i_txt = Textbox(fig; stored_string="1", validator=Int)
window_choose[1,4] = btn_next = Button(fig; label="+")
window_choose[1,5] = btn_next_on = Button(fig; label="⟩")
window_viz[1,1:3] = win_mids_ax = Axis(fig;
                                      aspect=DataAspect(),
                                      xticklabelsize=12, yticklabelsize=12)
mids_α_sl = Slider(fig; range=0:0.01:1)
mid_i_α_sl = Slider(fig; range=0:0.01:1)
img_α_sl = Slider(fig; range=0:0.01:1)
frame_i_sl = Slider(fig; range=1:1)
window_viz[2,1] = grid!( [Label(fig, "win. midlines") mids_α_sl;
                            Label(fig, "this midline") mid_i_α_sl;
                            Label(fig, "image")    img_α_sl;
                            Label(fig, "frame")    frame_i_sl] )
window_viz[2,2] = Label(fig, "lock")
window_viz[2,3] = win_ax_lock = Toggle(fig)
set_close_to!( mids_α_sl, 0.1 )
set_close_to!( mid_i_α_sl, 0.8 )
set_close_to!( img_α_sl, 0.8 )
hidedecorations!(win_mids_ax)

clamped_get(v,i) = v[clamp(i,firstindex(v),lastindex(v))]


# The last explicitly chosen window
window_i = async_latest(Observable(1))
# # Window currently displayed. My differ from `window_i on hover.
# window_i_disp = Observable(1)
onany(vis_data, window_i_txt.stored_string) do vd, str
    i = tryparse(Int, str)
    if window_i[] == i
        # Avoid loop with textbox updates.
    elseif i === nothing
        # Invalid text. Leave selected window as is.
    else
        clamped = clamp( i, 1, length(vd.windows ))
        window_i[] = clamped
    end
end
on(window_i) do i
    MakieLayout.set!(window_i_txt, string(i))
end
on(vis_data) do vd
    window_i[] = 1
end

window = @lift $vis_data.windows[$window_i]
on(window) do w
    frame_i_sl.range = w
    set_close_to!(frame_i_sl, frame_i_sl.value[])
end

win_mid_clines = lift(vis_data, window, mids_α_sl.value) do vd, w, α
    pts = replace(permutedims(vd.mids[w,begin:end]), missing=>Elegans.missingpoint)
    colors = repeat(w', size(pts,1))
    (; p = nan_separated_cols(pts), c = nan_separated_cols(colors),
       cmap = cgrad(:inferno, alpha=α) )
end

noframe_img = fill(coloralpha(colorant"black"), 151, 151)
noframe_imgdata = (; x = -75.0:75.0, y = -75.0:75.0, img = noframe_img)
img_data = Node(noframe_imgdata)

function set_image(vd, i, img, α)
    # Makie doesn't accept missings, or empty data (issue #533).
    # But `NaN`s work for plotting nothing.
    xi = coalesce(vd.traj.x[i], NaN)
    yi = coalesce(vd.traj.y[i], NaN)
    img_data[] = (; x = xi .+ (-75:75), y = yi .+ (-75:75),
                    img = alphacolor.(img, α) )
end

# set_no_image() = (img_data[] = noframe_imgdata)

waiting_for_img = Node(false)
onany(vis_data, frame_i_sl.value, img_α_sl.value) do vd, i, α
    function set_current_frame()
        img = get_frame(vd.vcache, i)
        # Transpose due to Makie bug #524: images by default plot like heatmaps,
        # with axis conventions that require manual transposition and y-axis flip to get
        # "image-like" axes
        set_image(vd, i, permutedims(img), α)
        #autolimits!(win_img_ax)
    end
    if Elegans.isframeready(vd.vcache, i)
        set_current_frame()
    else
        task = Task(set_current_frame)
        schedule(task)
        sleep(0.1)
        istaskdone(task) || (waiting_for_img[] = true)
        @async begin
            try
                wait(task)
            catch e
                showerror(stderr, e, catch_backtrace())
                println(stderr)
            finally
                waiting_for_img[] = false
            end
        end
    end
end

# trigger first window update to initialize image data and frame slider
window_i[] = 1
#MakieLayout.set!(window_i_txt, "1")

syncimage!( win_mids_ax, img_data )
text!(win_mids_ax, "Loading", position = @lift(($img_data.x[1], $img_data.y[1])),
    color = @lift $waiting_for_img ? :white : :transparent)
clines!( win_mids_ax, win_mid_clines)
lines!( win_mids_ax, @lift(replace($vis_data.mids[$(frame_i_sl.value),begin:end],
                                    missing=>Elegans.missingpoint)),
                    color = @lift(alphacolor(colorant"white", $(mid_i_α_sl.value))) )
onany(img_data, win_mid_clines) do _,_;
    win_ax_lock.active[] || autolimits!(win_mids_ax)
end


on_color = @lift RGBAf0(0,0,0,$(on_α.value))
off_color = @lift RGBAf0(0,0,0,$(off_α.value))
window_i_color = colorant"red"
window_clicked_color = window_i_color

selected_rect = select_rectangle(ax.scene)
# `active_rect` is either `selected_rect` or `nothing`
active_rect = Observable{Any}(selected_rect[])
on(selected_rect) do rect
    active_rect[] = rect
end

# indices of lines passing through selected_rect
inside_rect = @lift filter(i -> $active_rect === nothing || any(
                    (pt in $selected_rect) for pt in Point2.(s,$vis_data.y[:,i])),
                axes($vis_data.y,2))
# indices of samples from selected (on) and other (off) lines
on_i = Observable(Int[])
off_i = Observable(Int[])

function change_window(f)
    curr_i = window_i[]
    curr_i === nothing && return
    new_i = f( curr_i, on_i[] )
    window_i[] = new_i
    #MakieLayout.set!(window_i_txt, string(new_i))
end

on(btn_prev.clicks) do _; change_window((i,_)->max(i-1,1)) end
on(btn_next.clicks) do _; change_window((i,_)->min(i+1,length(vis_data[].windows))) end

on(btn_prev_on.clicks) do _
    change_window() do i, on_i
        prev_ind = searchsortedlast( on_i, i-1 )
        on_i[max(prev_ind,1)]
    end
end
on(btn_next_on.clicks) do _
    change_window() do i, on_i
        next_ind = searchsortedfirst( on_i, i+1 )
        on_i[min(next_ind,length(on_i))]
    end
end

function enable_button!(button, enabled = true, discolor = RGBf0(0.5,0.5,0.5))
    defaults = MakieLayout.default_attributes(Button, fig.scene).attributes
    for sym in (:buttoncolor, :buttoncolor_hover, :buttoncolor_active)
        button.:($sym) = enabled ? defaults[sym][] : discolor
    end
end

onany(window_i, on_i) do i, on
    @assert issorted(on)
    enable_button!(btn_prev_on, !isempty(on) && i > first(on))
    enable_button!(btn_next_on, !isempty(on) && i < last(on))
    enable_button!(btn_prev, i > 1)
    enable_button!(btn_next, i < length(vis_data[].windows))
end

ind2data(vd,i) = nan_separated_cols(s, @view vd.y[:,i])
on_data = lift(ind2data, vis_data, on_i)
off_data = lift(ind2data, vis_data, off_i)

on_off_e = onany(inside_rect, n_on, n_off, time_range_i) do inside_rect, n_on, n_off, time_range_i
    all_on = intersect( time_range_i, inside_rect )
    all_off = setdiff( time_range_i, all_on )

    on_i[] = sample(all_on, min(n_on, length(all_on)); replace=false) |> sort!
    off_i[] = sample(all_off, min(n_off, length(all_off)); replace=false) |> sort!
end

# async_latest causes trouble with Observables 0.3.3
# TODO: use it when GLMakie supports Observables 0.4
pdf_rect = #async_latest(
                latest_change(
                    lift( (active,rect) -> (active ? rect : nothing),
                          hist_sel_filter.active, active_rect, typ=Any )
                )
            #)

# reference pdfs are computed from `stage_time_span` rather than `time_range`
# to avoid double-update when changing stage/cam due to change of frame range.
refpdfs = lift(vis_data, latest_change(stage_time_span),
                #async_latest(hist_time_filter.active),
                hist_time_filter.active,
                pdf_rect
              ) do vd, tspan, tfilter, rect
    framerange = tfilter ? stage_timespan_to_framerange(tspan, vd.irange) : vd.irange
    pdfs = ref_pdf_estimates( allcams_data, s, vd.f, vd.cam_i, vd.stage_i,
                              framerange, rect;
                              edges=vd.y_edges )
    (; vd.y_edges, pdfs )
end

function pointwise_finmean(v)
    ax = only(unique(axes.(v)))
    [NaNMath.mean(filter(isfinite, [r[i] for r in v])) for i in CartesianIndices(ax)]
end

hmdata = lift(refpdfs, hm_type_menu.i_selected) do refpdfs, i
    pdfs, y_edges = refpdfs.pdfs, refpdfs.y_edges
    z = i == hm_types.this ?
            pdfs[cam_menu.i_selected[]] :
        i == hm_types.Δ ?
            pdfs[cam_menu.i_selected[]] .- mean(pdfs) :
        i == hm_types.logratio ?
            log2.(pdfs[cam_menu.i_selected[]] ./ mean(pdfs)) :
        i == hm_types.mean ?
            mean(pdfs) :
            #pointwise_finmean(pdfs) :
            error("i = $i")
    HMData(s_edges, y_edges, z)
end


pdf_crange = lift(hmdata, hm_type_menu.i_selected,
                  crange_toggle.active, crange_sl.value) do d, i, globalmax, v
    frac = 2.0^-v
    crange = if i == hm_types.logratio
        m = reduce(max, abs.(filter(isfinite, d.z)), init=0.0)
        (-m*frac, m*frac)
    else
        m = frac * (globalmax ? maximum(reduce(max, filter(isfinite,r), init=0.0) for r in refpdfs[].pdfs)
                              : reduce(max, filter(isfinite, d.z), init=0.0))
        if i ∈ (hm_types.this, hm_types.mean)
            (0.0, m)
        elseif i == hm_types.Δ
            (-m, m)
        else
            error("i = $i")
        end
    end
end

diverging = :diverging_bwr_40_95_c42_n256 # Alias missing from older `PlotUtils` versions
hm_cmap = lift(hm_type_menu.i_selected) do i
    i ∈ (hm_types.Δ, hm_types.logratio) ? cgrad(diverging) : cgrad(:viridis)
end


# trigger update
active_rect[] = nothing

hm = heatmap!(ax, hmdata, colormap=hm_cmap, colorrange = pdf_crange)

cb_layout[2,1] = cb = Colorbar(fig, hm, width=20)

on(hmdata) do _
    if hm_type_menu.i_selected[] ∈ (hm_types.Δ, hm_types.logratio)
        hm.colorrange = (-1,1) .* maximum(abs.(hm.colorrange[]))
    end
end


lines!(ax, on_data, color=on_color)
lines!(ax, off_data, color=off_color)
lines!(ax, s, @lift($vis_data.y[:,$window_i]), color=window_i_color, linewidth=2)

selected_rect_color = @lift RGBAf0(1,1,1, ($active_rect !== nothing) * 0.75)
deregister_interaction!(ax, :rectanglezoom)
poly!(ax, selected_rect, strokecolor = selected_rect_color, color=nothing, linestyle=:dot)
pt = select_point(ax.scene)

function enclosing_bin(edges, x)
    i = searchsortedfirst(edges,x)
    (from = i == firstindex(edges)  ? -Inf : edges[i-1],
     to   = i == lastindex(edges)+1 ?  Inf : edges[i])
end
function enclosing_2d_bin(edges_x, edges_y, pt)
    x_bin = enclosing_bin(edges_x, pt[1])
    y_bin = enclosing_bin(edges_y, pt[2])
    FRect( x_bin.from, y_bin.from, x_bin.to-x_bin.from, y_bin.to-y_bin.from )
end

allfinite(rect::Rect) = all(all(isfinite,p) for p in coordinates(rect))
pt_e = on(pt) do p
    # if point is near a rect vertex, assume a rectangle selection was
    # just performed and skip the point update.
    # TODO: find cleaner way to distinguish point from rect selection
    if all(!isapprox(p,q) for q in decompose(Point2{Float64}, selected_rect[]))
        bin = enclosing_2d_bin(s_edges, vis_data[].y_edges, p)
        allfinite(bin) && (selected_rect[] = bin)
    end
end

function index_from_mousepos( fig, upperbounds, dim=1 )
    pos = mouseposition(fig)[dim]
    i = searchsortedfirst(upperbounds, pos)
    min( i, lastindex(upperbounds) )
end

function hovered_index( axevents, fig, ub_obs, dim=1 )
    @assert issorted(ub_obs[])
    res = Observable(firstindex(ub_obs[]))
    handle = onmouseover(axevents) do evt
        res[] = index_from_mousepos( fig, ub_obs[], dim )
    end
    res, handle
end

axevents = addmouseevents!(ax.scene)
onmouserightup(axevents) do evt
    active_rect[] = nothing
end

spevents = addmouseevents!(speed_plot.scene)
onmouserightup(spevents) do evt
    reset_time_sel(vis_data[])
end

windows_upperbounds = @lift(last.($vis_data.windows))

hovered_bin, _ = hovered_index( axevents, ax.scene, Ref(s_edges[2:end]) )
hovered_window, _ = hovered_index( spevents, speed_plot.scene, windows_upperbounds )

hovering = Observable(false)

# Set the window to the hovered window while hovering, reset to last
# chosen window when exiting
explicit_window = Observable(1) # last window selected before hover, or by click
onmouseenter(spevents) do _
    hovering[] = true
    explicit_window[] = window_i[]
end
onmouseleftdown(spevents) do evt
    i = index_from_mousepos( speed_plot.scene, windows_upperbounds[] )
    window_i[] = explicit_window[] = i
end
on(hovered_window) do wi
    window_i[] = wi
end
onmouseout(spevents) do _
    window_i[] = explicit_window[]
    hovering[] = false
end

vlines!(ax, @lift(s[$hovered_bin]), color=RGBAf0(0,0,0,0.5), linestyle=:dot)
vlines!(speed_plot, @lift(middle.($vis_data.windows[$inside_rect])), color=RGBAf0(0,0,0,0.05))
vlines!(speed_plot, @lift(middle.($vis_data.windows[$on_i])), color=RGBAf0(0,0,0,0.5))
vlines!(speed_plot, @lift(middle(clamped_get($vis_data.windows, $window_i))),
                    color = window_i_color,
                    #linestyle = @lift(ifelse($hovering, :dot, :solid)),
                    linewidth = 3 )
vlines!(speed_plot, @lift(middle(clamped_get($vis_data.windows, $explicit_window))),
                    color = @lift(coloralpha(window_i_color, $hovering)),
                    linestyle=:dot, linewidth=3)
# This needs to be after all plot commands to `speed_plot` to ensure it happens
# after all plot data updates
on(vis_data) do vd
    # time_sel[] = Point2.( [first(vd.irange),last(vd.irange)], speed_lims[] )
    autolimits!(speed_plot)
end

lines!(hist1d_ax, @lift(Point2.($hmdata.z[$hovered_bin,:], midpoints($hmdata.y))))
#lines!(hist1d_ax, @lift($z[$hovered_bin,:]), lift(midpoints,edges))
_extrema(x, default) = isempty(x) ? default : extrema(x)
set_hist1d_lims(d) = limits!(hist1d_ax,
                            reverse(_extrema(filter(isfinite,d.z), (0,1e-12))),
                            extrema(d.y))
set_hist1d_lims(hmdata[])
on(set_hist1d_lims, hmdata)
hist1d_ax.ylabel = "log₁₀(GV⋅len⁻²)"

display(fig)
