function out = combine_points(points,values,d)
% ������ļ��С����ֵd����ϲ���������
% ע��points��ʾ������������꣬d��ʾ������ֵ

points_num = length(points(:,1));
flag = zeros(points_num,1);
out = zeros(points_num,2);
out_index = 0;
d = d^2;
while true
    list = find(flag == 0);
    if isempty(list)
        break;
    else % ����length(list)=1���������ֻʣ�����һ����
        distance = sum((points(list,:) - repmat(points(list(1),:),length(list),1)).^2,2);
        near_index = list(distance < d);
        if length(near_index) == 1
            flag(near_index) = 1;
            out_index = out_index + 1;
            out(out_index,:) = points(near_index,:);
        else
            flag(near_index) = 1;
            % ֻ���ǵ�ļ���λ��
            % points(near_index(1),:) = mean(points(near_index,:));
            % ֻ�������ֵ�ĵ�
            maxvalue_index = find(values(near_index)==max(values(near_index)),1);
            points(near_index(1),:) = points(list(maxvalue_index(1)),:);
            
            % ���µ�ʱ�������⣬��ɨ����ͨ��֧�����㲻ͬ��ͨ��֧�ľ���
            flag(near_index(1)) = 0;%���µ㣬���²����ϲ�����
        end
    end
end
out = out(1:out_index,:);
end