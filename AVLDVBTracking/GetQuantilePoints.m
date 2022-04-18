function quantile = GetQuantilePoints(Folder, quantile_index) 
CenterlineFolder = [Folder 'Wave\centerline_GCaMP\'];
CenterlineFiles = dir([CenterlineFolder '*.mat']);
quantile = zeros(length(CenterlineFiles),2);
for i = 1:length(CenterlineFiles)
    disp(['Processing ' num2str(i) '/' num2str(length(CenterlineFiles))]);
    Centerline = load([CenterlineFolder num2str(i) '.mat']);
    quantile(i,1) = Centerline.centerline(quantile_index,2);
    quantile(i,2) = Centerline.centerline(quantile_index,1); 
end
output_name = [Folder sprintf('neuron_pos\green\centerline %02d',quantile_index) '.txt'];
fid = fopen(output_name,'w');  
for i = 1:length(CenterlineFiles)
    fprintf(fid,'%d    %d\n',quantile(i,1),quantile(i,2));
end
    