function ends_alignment_mask(ends)
    out = Vector{Bool}(undef,length(ends))
    out[1] = false   # Just for consistency. Not necessary for correct output.
    prev = ends[1]
    for i = 2:length(ends)
        curr = ends[i]
        curr_reversed = curr[[2,1]]
        out[i] = out[i-1] ⊻ (mean(norm.(curr .- prev)) > mean(norm.(curr_reversed .- prev)))
        prev = curr
    end
    out
end

function align_ends!(lines)
    ends = [(l[1], l[end]) for l in lines]
    switch_ends = ends_alignment_mask(ends)
    for i in eachindex(ends)
        switch_ends[i] && (lines[i] = reverse(lines[i]))
    end
    lines
end
