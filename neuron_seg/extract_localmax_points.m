function localmax_points = extract_localmax_points(localmax_img, worm_img,conncomp_region)
% 提取局部最大值点

start_r = conncomp_region(1);
start_c = conncomp_region(2);
height = conncomp_region(3);
width = conncomp_region(4);

localmax_region = zeros(size(localmax_img));
localmax_region(start_r:start_r+height-1,start_c:start_c+width-1)...
    = localmax_img(start_r:start_r+height-1,start_c:start_c+width-1);

localmax_conncomps = bwconncomp(localmax_region,8);
for i=1:localmax_conncomps.NumObjects
    pixelslist = localmax_conncomps.PixelIdxList{i};
    localmax_values = worm_img(pixelslist);
    % maxvalue_index = find(localmax_values == max(localmax_values),1);
    maxvalue_index = find(localmax_values > 0.85*max(localmax_values));
    localmax_region(pixelslist) = 0;
    localmax_region(pixelslist(maxvalue_index)) = 1;
end
localmax_points = find(localmax_region==1);
end