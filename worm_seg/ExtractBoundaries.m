function [left,right,centerline] = ExtractBoundaries(bw_img,centerline)
% ���߳��ֵͼ��bw_img����ȡ�߳����ߣ����ұ�Ե�������ߣ�
% left, right, center�ֱ��ʾ�߳����Ե���ұ�Ե��������
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

% % % make worm thinner (ʹͼ���Ե���ӽ������Ե��
% se = strel('disk',2);
% bw_img = imerode(bw_img, se);

% ���������߼������ұ�Ե
% width = 50;
% centerline = centerline(2:end-1,:); % ȥ�����Ե�ĵ㣬ȷ�������׼ȷ��
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
    dir = centerline(i+1,:) - centerline(i,:);% ��β����ǰ����
    dir = dir/norm(dir);
    dir_p = [-dir(2) dir(1)];% �����߷���
    %�ط����������������������������������ұ�Ե(dorsal)�������������Ե��ventral)
    %�ұ�Ե��
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
    
    %���Ե��
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