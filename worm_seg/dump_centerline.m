function dump_centerline(centerline,filename)
% save centerline with filename as txt format

Interval = 2;
point_index = 1:size(centerline,1);

downsample_centerline = centerline(1:Interval:end,:);
downsample_point_index = point_index(1:Interval:end);

% add tail into downsample centerline
if downsample_point_index(end) < size(centerline,1)
    downsample_centerline = [downsample_centerline; centerline(end,:)];
end

% write downsample centerline into file
fid = fopen(filename,'w+t');
for t = 1:size(downsample_centerline,1)
    fprintf(fid,'%d    %d\n',downsample_centerline(t,2),downsample_centerline(t,1));
end
fclose(fid);
end

