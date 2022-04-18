function [left,right,centerline] = ExtractBoundaries(bw_img,centerline)
% 从线虫二值图像bw_img中提取线虫曲线（左右边缘，中心线）
% left, right, center分别表示线虫左边缘，右边缘和中心线
% close all;

%  image processing to fill holes
% se = strel('disk', 3); % fill small holes inside worm region
% bw_img = (bw_img > 0);
% bw_img = imclose(bw_img,se);

margin = 5; % clear border of image
[height, width] = size(bw_img);
bw_img(1:margin, :) = 0;
bw_img((height-margin+1):height, :) = 0;
bw_img(:,1:margin) = 0;
bw_img(:,(width-margin+1):width) = 0; 

% % % make worm thinner (使图像边缘更接近肌肉边缘）
% se = strel('disk',2);
% bw_img = imerode(bw_img, se);

% 根据中心线计算左右边缘
% width = 50;
% centerline = centerline(2:end-1,:); % 去掉最边缘的点，确保计算的准确性
width = 1000;
PointsNum = length(centerline);
tail = centerline(PointsNum,:);
head = centerline(1,:);
boundary_image = inner_edge(bw_img);

left = zeros(PointsNum,2);
right = zeros(PointsNum,2);
left(PointsNum,:) = tail;
right(PointsNum,:) = tail;
left(1,:) = head;
right(1,:) = head;

for i=PointsNum-1:-1:2
    dir = centerline(i+1,:) - centerline(i,:);% 从尾部向前搜索
    dir = dir/norm(dir);
    dir_p = [-dir(2) dir(1)];% 中心线法线
    %沿法线向两边延伸获得左右轮廓，正方向是右边缘(dorsal)，负方向是左边缘（ventral)
    %右边缘点
    delta_x = 2*(dir(1)>0)-1;
    delta_y = 1-2*(dir(2)>0);
    for j=0:width        
        point = j*dir_p + centerline(i,:);
        p = int32(point);
        if boundary_image(p(1),p(2))==1 || boundary_image(p(1)+delta_x,p(2))==1 ||...
           boundary_image(p(1),p(2)+delta_y)==1 || boundary_image(p(1)+delta_x,p(2)+delta_y)==1 ||...
           boundary_image(p(1)-delta_x,p(2))==1 || boundary_image(p(1)-delta_x,p(2)+delta_y)==1
            right(i,:) = point + [0.5,0.5].*[delta_x,delta_y];
            break;
        end       
    end
    
    %左边缘点
    delta_x = -delta_y;
    delta_y = delta_x;
    for j=0:width
        point = j*(-dir_p) + centerline(i,:);
        p = int32(point);
        if boundary_image(p(1),p(2))==1 || boundary_image(p(1)+delta_x,p(2))==1 ||...
           boundary_image(p(1),p(2)+delta_y)==1 || boundary_image(p(1)+delta_x,p(2)+delta_y)==1 ||...
           boundary_image(p(1)-delta_x,p(2))==1 || boundary_image(p(1)-delta_x,p(2)+delta_y)==1
            left(i,:) = point + [0.5,0.5].*[delta_x,delta_y];
            break;
        end
    end
end

% local search to find the strong boundary

% smooth the left and right boundaries
% Interval = 1;
% P_Num = 50;
% left_part = left(1:Interval:PointsNum,:);
% right_part = right(1:Interval:PointsNum,:);
% left = spline_fitting_partition(left_part,P_Num);
% right = spline_fitting_partition(right_part,P_Num);

% based on centerline, to determine points in left and right boundaries


% % show boundarries
% imagesc(bw_img);axis image;colormap(gray);hold on;
% plot(left(:,2),left(:,1),'r.');
% plot(right(:,2),right(:,1),'g.');
% plot(centerline(:,2),centerline(:,1),'b-');

end