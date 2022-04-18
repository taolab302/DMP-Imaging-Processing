function in_list = points_in_list(points,list)
% �жϸ����㼯points�Ƿ���list�У�������
% ע��points��list����ʾ���index

if numel(points) == 0 || numel(list) == 0
    in_list = [];
else
    in_list = zeros(1,length(points));
    index = 0;
    for i=1:length(points)
        if ~isempty(find(list == points(i),1))
            index = index+1;
            in_list(index) = points(i);
        end
    end
    if index ~= 0
        in_list = in_list(1:index)';
    else
        in_list = [];
    end
end
end