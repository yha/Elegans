using LRUCache

struct VideoCache{F}
    path_f::F
    boundaries::Vector{Int}
    cache::LRU{Int,Any}

    function VideoCache( cam, datadir = Elegans.datadir, maxsize = 20 )
        path_f = Elegans.videopath_f(cam, datadir)
        _, boundaries = import_coords(cam, datadir)
        cache = LRU{Int,Any}( maxsize = maxsize )
        new{typeof(path_f)}( path_f, boundaries, cache )
    end
end

function get_frame(cache::VideoCache, i)
    fileno = searchsortedfirst(cache.boundaries, i) - 1
    file_offset = fileno == 0 ? 0 : cache.boundaries[fileno]
    # @show (fileno,i-file_offset)
    frames = get!(cache.cache, fileno) do
        println("Reading file $fileno")
        Elegans.read_video( cache.path_f, fileno )
    end
    frame = frames[i-file_offset]
end

nframes(cache::VideoCache) = cache.boundaries[end]
