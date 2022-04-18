function res = LoadCenterlineResults(filename)

points_num = 101;
fid = fopen(filename,'rb');
res = struct('length_error',false,'last_width',0,...
    'last_mean_length',0,'current_backbone',zeros(points_num,2),...
    'last_backbone',zeros(points_num,2),'actual_backbone',zeros(points_num,2));

res.length_error = fread(fid,1,'bool');
res.last_mean_length = fread(fid,1,'double');
res.last_width = fread(fid,1,'double');

actual_backbone_num = fread(fid,1,'int');
backbone = fread(fid,points_num*2,'double');
res.actual_backbone(:,1) = backbone(1:2:end);
res.actual_backbone(:,2) = backbone(2:2:end);

current_backbone_num = fread(fid,1,'int');
backbone = fread(fid,points_num*2,'double');
res.current_backbone(:,1) = backbone(1:2:end);
res.current_backbone(:,2) = backbone(2:2:end);

last_backbone_num = fread(fid,1,'int');
backbone = fread(fid,points_num*2,'double');
if ~isempty(backbone)
    res.last_backbone(:,1) = backbone(1:2:end);
    res.last_backbone(:,2) = backbone(2:2:end);
end

fclose(fid);
end