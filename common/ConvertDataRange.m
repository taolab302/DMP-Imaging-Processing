function range = ConvertDataRange(Sync_Struc,channel,data_range)
% Convert data range in channel (green/red) to range in (red/green)

channel = lower(channel);
sync_index = Sync_Struc.match_index;
if strcmp(channel,'g') == 1 || strcmp(channel,'green') == 1   %convert green to red
    range = sync_index(data_range);   
elseif strcmp(channel,'r') == 1 || strcmp(channel,'red') == 1 %convert red to green
    range = [1,1];
    range(1) = find(sync_index == data_range(1),1);
    range(2) = find(sync_index == data_range(2),1,'last');
    
    % if there is no corresponding green channel image for the last red channel image
    if isempty(range(2))
        range(2) = find(sync_index == data_range(2)-1,1,'last');
    end
end
end