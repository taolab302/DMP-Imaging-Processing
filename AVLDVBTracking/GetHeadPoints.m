function head = GetHeadPoints(Folder) 
WormRegionPosFile = load([Folder 'Wave\WormRegionPos.mat']);
head = zeros(length(WormRegionPosFile.centerline_start),2);
for i = 1:length(head)
    disp(['Processing ' num2str(i) '/' num2str(length(head))]);
    head(i,1) = WormRegionPosFile.centerline_start(i,2);
    head(i,2) = WormRegionPosFile.centerline_start(i,1);
end
output_name = [Folder 'neuron_pos\green\centerline 01.txt'];
fid = fopen(output_name,'w');  
for i = 1:length(head)
    fprintf(fid,'%d    %d\n',head(i,1),head(i,2));
end
    