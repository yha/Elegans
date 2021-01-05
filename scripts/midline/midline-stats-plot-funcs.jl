using Statistics, StatsBase, KernelDensity
using Plots, LaTeXStrings
using NaNMath


function allvars_hm( vars, t, windows, isroam, summ_str )
    midvars = vars[:,size(vars,2)÷2+1]

    hm = heatmap(t, first.(windows), log2.(vars./midvars), yaxis=false, title=summ_str)
    plot!(  hm, fill(0.5,size(windows)), first.(windows), label="",
            lz=log10.(midvars), colorbar_entry=false, c=:grays, lw=2 )
    roamwins = filter(win->any(coalesce.(isroam[win],false)), windows)
    isempty(roamwins) || hline!(first.(roamwins), ls=:dot, lc="black", lw=2, label="")
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


function make_all_gv_figures( gv, t, windows, edges, summ_str, isroam; sample=axes(gv), ftype="png", path="." )
    cd(path) do
        savefig(allvars_hm(gv,t,windows,isroam,summ_str), "midpoint GV $summ_str heatmap.$ftype")
        savefig(allvars_plot(gv,t,summ_str,sample), "midpoint GV $summ_str sample plot.$ftype")
        savefig(vars_hist_hm(gv,t,edges,summ_str), "midpoint GV $summ_str hist.$ftype")
        savefig(vars_kde_hm(gv,t,edges,summ_str), "midpoint GV $summ_str kde.$ftype")
    end
end
