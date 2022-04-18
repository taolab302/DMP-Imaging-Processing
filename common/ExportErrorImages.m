function ExportErrorImages(Image_Folder,file_op)
% Export incorrectly calculated image data

Error_Folder = [Image_Folder 'error'];
if ~exist(Error_Folder,'dir')
    mkdir(Error_Folder);
end
Error_Image_Folder = [Error_Folder '\Image'];
if ~exist(Error_Image_Folder,'dir')
    mkdir(Error_Image_Folder);
end
Error_Image_Folder = [Error_Image_Folder '\'];

offset = 1;
error_list = load([Image_Folder 'list.txt']);
error_list = error_list + offset - 1;
error_list = unique(error_list);
if isempty(error_list)
    return;
end

% copy incorrectly calculated image data into error folder
Image_Folder_Prefix = 'Image\';
% Image_Folder_Prefix = 'RFP_Map\';
image_names = dir([Image_Folder Image_Folder_Prefix '*.tiff']);
for i=1:length(error_list)
    image_index = error_list(i);
    disp(['Copy File: ' num2str(image_index)]);
    imagename = [Image_Folder Image_Folder_Prefix char(image_names(image_index).name)];
    copyfile(imagename,[Error_Image_Folder char(image_names(image_index).name)]);
end

% empty the error list
if strcmp(file_op,'clr') ==  1
    file = fopen([Image_Folder 'list.txt'],'wt');
    fprintf(file,'\n');
    fclose(file);
end
end

