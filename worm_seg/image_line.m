function output = image_line(im,points)
% 在图像im中连接p1和p2两点(x,y)

p1 = points(1,:);
p2 = points(2,:);
[height,width] = size(im);
output = im;

if (p1(1)<0 || p1(1)>width) || (p1(2)<0 || p1(2)>height)||...
   (p2(1)<0 || p2(1)>width) || (p2(2)<0 || p2(2)>height)
    disp('More than one points out of range');
    return;
end

% 连接线段[p1,p2]
point1 = round(p1);
point2 = round(p2);
if (point1(1) == point2(1))
    min_y = min(point1(2),point2(2));
    max_y = max(point1(2),point2(2));
    output(min_y:max_y,point1(1)) = 1;
else
    k = (point2(2)-point1(2))/(point2(1)-point1(1));
    output(point1(2),point1(1)) = 1;
    output(point2(2),point2(1)) = 1;
    if abs(k)<=1
        delta_x = 2*(point2(1)>point1(1))-1;
        k = k/delta_x; %修正方向
        y = point1(2);
        for x=point1(1):delta_x:point2(1)-delta_x
            output(round(y),x) = 1;
            y = y+k;
        end
    else
        delta_y = 2*(point2(2)>point1(2))-1;
        k = k/delta_y; %修正方向
        x = point1(1);
        for y=point1(2):delta_y:point2(2)-delta_y
            output(y,round(x)) = 1;
            x = x+1/k;
        end
    end
end
end