function [seed_img,new_seed_points] = create_seeds(seed_points,img_region)
% 根据图像区域，生成种子图像

start_row = img_region(1);
start_col = img_region(2);
height = img_region(3);
width = img_region(4);
margin = img_region(5);

seed_img = zeros(height+2*margin,width+2*margin);
new_seed_points = zeros(size(seed_points));
new_seed_points(:,1) = seed_points(:,1)-start_row+margin+1;
new_seed_points(:,2) = seed_points(:,2)-start_col+margin+1;

for i=1:length(seed_points(:,1))
    seed_img(new_seed_points(i,1),new_seed_points(i,2)) = 1;
end
end