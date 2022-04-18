function edgeList = EdgeList(data)
% 本程序用来求取二值图的边缘链表
% data为二值图像，1为目标，0为背景
% edgeList 为data的边缘链表2*N，其首链表和末链表为同一边缘点

direction = [0 1; -1 1; -1 0; -1 -1; 0 -1; 1 -1; 1 0; 1 1];
direction = [direction; direction];

% 扩展数据，使其周边像素为0
im = zeros(size(data,1) + 2, size(data,2) + 2);
im(2:(end-1), 2:(end-1)) = data;

%查找起始点
[x, y] = find(im == 1);
xori = x(1); yori = y(1);

edgex = xori;
edgey = yori;

% 设置起始查找方向
i = 6;
while 1
    if (1 == im(edgex(end) + direction(i, 1), edgey(end) + direction(i, 2)))
        edgex = [edgex, edgex(end) + direction(i, 1)];
        edgey = [edgey, edgey(end) + direction(i, 2)];

        offset = 2 - mod(i, 2);
        i = mod( (i - offset) - 1, 8) + 1 ;     
        if(edgex(1) == edgex(end) && edgey(1) == edgey(end))
            break;
        end
    else
        i = i+1;
    end    
end

edgeList = [edgex-1; edgey-1];
end
