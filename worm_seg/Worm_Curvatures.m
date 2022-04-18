function curvatures = Worm_Curvatures(Centerline_Folder,frame_seq,SkipList)
% calculate worm curvature

centerlines = dir([Centerline_Folder, '*.mat']);


if strcmp(frame_seq,'all')
    Num = length(centerlines);
    frame_seq = zeros(1,Num);
    for i = 1:Num
        frame_seq(i) = str2double(centerlines(i).name(1:end-4));
    end
    frame_seq = sort(frame_seq);
else
    Num = length(frame_seq);
end
    
for i = 1:Num
    if ~isempty(find(SkipList == (i-1), 1))
        centerline_num = size(curvatures,2); % assume the first iisn't in skiplist
        curvatures(i,:) = nan(1,centerline_num);
        continue;
    end
    
	centerline_data = load([Centerline_Folder num2str(frame_seq(i)) '.mat']);
    centerline = centerline_data.centerline;
	if i==1
		curvatures = zeros(Num,length(centerline)); % allocate spaces
	end
	curvatures(i,:) = Compute_Curvature(centerline);
end
% % save worm regions and positions
% save('WormCurvature.mat','curvatures');
end