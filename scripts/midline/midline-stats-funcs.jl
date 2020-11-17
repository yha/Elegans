# using Revise
# using DrWatson
# @quickactivate
# ##
using Elegans
using Images, ImageFiltering
using ImageFiltering: centered  # conflict with GeometryBasics
using Statistics, StatsBase, KernelDensity
using StaticArrays, OffsetArrays
using LinearAlgebra
using GeometryBasics
using Plots, LaTeXStrings
using Missings, NaNMath


##

function load_cam(root, ex, cam, contours_path, midpoints_path, contour_method=Thresholding(1.0,0.34))
    relcam = joinpath(ex, cam)
    campath = joinpath(root, relcam)
    @assert isdir(campath)

    contours, contours_file, vcache = init_contours(
                            relcam, root, contour_method, contours_path)

    traj = import_and_calc(relcam, 3, root)

    mids, midfile = Elegans.init_midpoints(ex, cam, traj, contours, midpoints_path; contour_method)

    (;traj, contours, mids, vcache, contours_file, midfile)
end


# used instead of skipmissing as a workaround for Statistics.jl issue #50
nomiss(v) = disallowmissing(filter(!ismissing,v))
_cov(x::AbstractVector{P}) where {N,T,P<:Point{N,T}} = isempty(x) ? fill(T(NaN),SMatrix{2,2}) : cov(x)

function framerange( ex, cam, stage_i; stages=loadstages() )
    stage_boundaries = stages[ex][cam]
    return stage_boundaries[stage_i]+1:stage_boundaries[stage_i+1]
end

function midline_covs( ex, cam, mids, stage_i; winlen=60, dwin=20, stages=loadstages() )

    irange = framerange( ex, cam, stage_i; stages )
    @info "$cam, stage $stage_i: frames $irange"
    midpts = mids(irange)

    t = mids.t
    windows = [irange[i0+1:i0+winlen] for i0 in 0:dwin:length(irange)-winlen]
    @time covs = [_cov(nomiss(midpts[i,j])) for i in windows, j in eachindex(t)]

    (;covs, windows, t, irange)
end

const newaxis = [CartesianIndex()]
const newaxis_c = OffsetVector([CartesianIndex()],0:0)
function normal_speeds(midpts)
    d = centered([-1/2,0,1/2]) # derivative

    midpts_nan = replace(midpts, missing=>Elegans.missingpoint)
    # tangents pointing towards tail
    tangents = normalize.(imfilter(midpts_nan, d[newaxis_c,:]))
    # normals pointing to the worm's right
    nrm = [[-p[2],p[1]] for p in tangents]

    v = imfilter(midpts_nan, d[:,newaxis_c])
    #v_t = dot.(tangents, v)
    v_n = dot.(nrm, v)
end



##


function allvars_hm( vars, t, windows, isroam, summ_str )
    midvars = vars[:,size(vars,2)÷2+1]

    hm = heatmap(t, first.(windows), log2.(vars./midvars), yaxis=false, title=summ_str)
    #plot!(hm, [0.5, 0.5], [3.5e5, 4e5], lz=[100,200], colorbar_entry=false)
    plot!(  hm, fill(0.5,size(windows)), first.(windows), label="",
            lz=log10.(midvars), colorbar_entry=false, c=:grays, lw=2 )
    hline!(filter(i->coalesce(isroam[i],false),
                    first.(windows)), ls=:dot, lc="black", lw=2, label="")
    plot( #plot(midvars, first.(windows), xscale=:log10, legend=false, xflip=true),
          plot(log10.(midvars), first.(windows), legend=false, xflip=true),
          hm,
          layout=grid(1,2,widths=(0.2,0.8)), link=:y )
end

function allvars_plot( vars, t, summ_str, sample=axes(vars) )
    plot(t, view(vars,sample...)', c="black", title="variances \n $summ_str",
            la=0.025, label=false, ylabel="GV", yscale=:log10, legend=:bottomright)
    plot!( t, mapslices(NaNMath.mean, vars, dims=1)[1,:], lw=2, label="mean" )
    plot!( t, exp.(mapslices(NaNMath.mean, log.(vars), dims=1)[1,:]), lw=2, label="geom. mean" )
    plot!( t, mapslices(NaNMath.median, vars, dims=1)[1,:], c=3, lw=2, label="median" )
    #savefig("midpoint-vars-all $summ_str.png"); current()
end

function vars_hist_hm( vars, t, edges, summ_str )
    hists = [fit(Histogram, log10.(filter(isfinite,col)), edges) for col in eachcol(vars)]

    heatmap(t, StatsBase.midpoints(edges), reduce(hcat, h.weights for h in hists), legend=false, cbar=true,
            xlabel=L"s", ylabel=L"\log_{10}(\mathrm{GV})",
            title="variance histograms \n $summ_str")
    plot!( t, log10.(mapslices(NaNMath.mean, vars, dims=1)[1,:]), c=1, lc="white", ls=:dash )
    plot!( t, mapslices(NaNMath.mean, log10.(vars), dims=1)[1,:], c=1, lc="white", ls=:dash )
    plot!( t, mapslices(NaNMath.median, log10.(vars), dims=1)[1,:], c="lightgray", ls=:dash )
end

function vars_kde_hm( vars, t, edges, summ_str )
    k = [kde(filter(isfinite,log10.(col))) for col in eachcol(vars)]
    ik = InterpKDE.(k)

    heatmap(t, edges, [pdf(i,x) for x in edges, i in ik], legend=:bottomleft,
            xlabel=L"s", ylabel=L"\log_{10}(\mathrm{GV})",
            title="variance KDE \n $summ_str")
    plot!( t, mapslices(log10∘NaNMath.mean, vars, dims=1)[1,:], c="gray", label="mean" )
    plot!( t, mapslices(NaNMath.mean, log10.(vars), dims=1)[1,:], c="lightgray", label="g.m." )
    plot!( t, mapslices(NaNMath.median, log10.(vars), dims=1)[1,:], c=1, label="median" )
end

##

function make_all_gv_figures( gv, t, windows, edges, summ_str, isroam; sample=axes(gv), ftype="png", path="." )
    cd(path) do
        savefig(allvars_hm(gv,t,windows,isroam,summ_str), "midpoint GV $summ_str heatmap.$ftype")
        savefig(allvars_plot(gv,t,summ_str,sample), "midpoint GV $summ_str sample plot.$ftype")
        savefig(vars_hist_hm(gv,t,edges,summ_str), "midpoint GV $summ_str hist.$ftype")
        savefig(vars_kde_hm(gv,t,edges,summ_str), "midpoint GV $summ_str kde.$ftype")
    end
end
