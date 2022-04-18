function out = combine_points(points,values,d)
% 若两点的间隔小于阈值d，则合并这两个点
% 注：points表示待操作点的坐标，d表示距离阈值

points_num = length(points(:,1));
flag = zeros(points_num,1);
out = zeros(points_num,2);
out_index = 0;
d = d^2;
while true
    list = find(flag == 0);
    if isempty(list)
        break;
    else % 包括length(list)=1的情况，即只剩下最后一个点
        distance = sum((points(list,:) - repmat(points(list(1),:),length(list),1)).^2,2);
        near_index = list(distance < d);
        if length(near_index) == 1
            flag(near_index) = 1;
            out_index = out_index + 1;
            out(out_index,:) = points(near_index,:);
        else
            flag(near_index) = 1;
            % 只考虑点的几何位置
            % points(near_index(1),:) = mean(points(near_index,:));
            % 只考虑最大值的点
            maxvalue_index = find(values(near_index)==max(values(near_index)),1);
            points(near_index(1),:) = points(list(maxvalue_index(1)),:);
            
            % 更新点时存在问题，可扫描连通分支，计算不同连通分支的距离
            flag(near_index(1)) = 0;%更新点，重新参与点合并操作
        end
    end
end
out = out(1:out_index,:);
end