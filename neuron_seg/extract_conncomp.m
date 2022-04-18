function [conncomp_img,points,region] = extract_conncomp(img,conncomp)
% 从图像img中提取连通分支conncomp对应的图像

[height,width] = size(img);
points_num = length(conncomp);
points = zeros(points_num,2);
points(:,2) = floor(conncomp/height)+1;     %列坐标
points(:,1) = conncomp - (points(:,2)-1)*height;%行坐标

margin = 11;
min_r = min(points(:,1));
min_c = min(points(:,2));
conncomp_height = max(points(:,1))-min_r+1;
conncomp_width  = max(points(:,2))-min_c+1;
region = [min_r,min_c,conncomp_height,conncomp_width,margin];

conncomp_img = zeros(conncomp_height+2*margin, conncomp_width+2*margin);
for i=1:points_num
    conncomp_img(points(i,1)-min_r+margin+1,points(i,2)-min_c+margin+1) =...
        img(points(i,1),points(i,2));
end
end