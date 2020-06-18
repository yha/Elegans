using ImageFiltering
using Images.KernelFactors: gaussian

function log_speed_ratios(end1, end2, ends_found, t, s)
    v1 = diff(imfilter(end1, gaussian(t), NA()))
    v2 = diff(imfilter(end2, gaussian(t), NA()))
    n1 = norm.(v1)
    n2 = norm.(v2)
    n1 = replace(spread([1; n1], ends_found), missing => 1)
    n2 = replace(spread([1; n2], ends_found), missing => 1)
    imfilter(log.(n1 ./ n2), gaussian(s), NA())
end
