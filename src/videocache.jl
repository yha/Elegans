using LRUCache

struct VideoCache
    well::Well
    boundaries::Vector{Int}
    cache::LRU{Int,Any}

    @deprecate VideoCache(well,datadir) VideoCache(Well(datadir, splitdir(well)...))
    function VideoCache( well::Well; traj = load_coords_and_size(well), maxsize = 20 )
        boundaries = boundaries_from_traj(traj)
        cache = LRU{Int,Any}( maxsize = maxsize )
        new( well, boundaries, cache )
    end
end

function globalframe_to_fileframe(cache::VideoCache, i)
    fileno = searchsortedfirst(cache.boundaries, i) - 1
    file_offset = fileno == 0 ? 0 : cache.boundaries[fileno]
    fileno, i - file_offset
end

function _fetch_video(cache::VideoCache, fileno)
    frames = get!(cache.cache, fileno) do
        @info "Reading video $fileno of $(cache.well)" _group=:videoread
        Elegans.load_video( cache.well, fileno )
    end
end

function get_frame(cache::VideoCache, i; unpad = true)
    fileno, frameno = globalframe_to_fileframe(cache, i)
    frames = _fetch_video(cache, fileno)
    frame = frames[frameno]
    if unpad
        # Frame should have odd size, but may be padded to even size
        # due to encoding limitations. Undo padding.
        frame = frame[1:end-mod((size(frame,1)+1),2),
                    1:end-mod((size(frame,2)+1),2)]
    end
    frame
end

function prefetch_framerange(cache::VideoCache, from, to)
    from_file, _ = globalframe_to_fileframe(cache, from)
    to_file, _   = globalframe_to_fileframe(cache, to)
    @info "Prefetching video files $(from_file)—$to_file"
    for k in from_file:to_file
        _fetch_video(cache, k)
        yield()
        #sleep(0.05)
    end
end 
prefetch_framerange(cache::VideoCache, r::AbstractUnitRange) = prefetch_framerange(cache, first(r), last(r))

nframes(cache::VideoCache) = cache.boundaries[end]

isframeready(cache::VideoCache, frame_index) = video_index(cache, frame_index) ∈ cache.cache.keyset

video_index(cache::VideoCache, frame_index) = _video_index(cache.boundaries, frame_index)
