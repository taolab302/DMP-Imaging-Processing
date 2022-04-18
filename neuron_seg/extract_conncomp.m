function [conncomp_img,points,region] = extract_conncomp(img,conncomp)
% ��ͼ��img����ȡ��ͨ��֧conncomp��Ӧ��ͼ��

[height,width] = size(img);
points_num = length(conncomp);
points = zeros(points_num,2);
points(:,2) = floor(conncomp/height)+1;     %������
points(:,1) = conncomp - (points(:,2)-1)*height;%������

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