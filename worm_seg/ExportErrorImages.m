function ExportErrorImages(Image_Folder)
% Export incorrectly calculated image data

Error_Folder = [Image_Folder 'error\'];
if ~exist(Error_Folder,'dir')
    mkdir(Error_Folder);
end

offset = 1;
error_list = load([Image_Folder 'list.txt']);
error_list = error_list + offset - 1;
error_list = unique(error_list);
if isempty(error_list)
    return;
end

worm_regions = load([Image_Folder 'WormRegionPos.mat']);
worm_regions = worm_regions.worm_regions;

% copy incorrectly calculated image data into error folder
Image_Folder_Prefix = 'worm_region\';
% Image_Folder_Prefix = 'RFP_Map\';
image_names = dir([Image_Folder Image_Folder_Prefix '*.tiff']);
for i=1:length(error_list)
    image_index = error_list(i);
    disp(['Copy File: ' num2str(image_index-1)]);
    imagename = [Image_Folder Image_Folder_Prefix num2str(image_index-1) '.tiff'];
    err_image = imread(imagename);

    err_image = err_image(worm_regions(image_index,1):worm_regions(image_index,2),worm_regions(image_index,3):worm_regions(image_index,4));
    imwrite(err_image,[Error_Folder  num2str(image_index-1) '.tiff']);
end

% % empty the error list
% if strcmp(file_op,'clr') ==  1
%     file = fopen([Image_Folder 'list.txt'],'wt');
%     fprintf(file,'\n');
%     fclose(file);
% end
end

