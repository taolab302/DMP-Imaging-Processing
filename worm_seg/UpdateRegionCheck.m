function UpdateRegionCheck(Image_Folder,wave_index)
%     map_imgFolder = [Image_Folder(1:end-1) '-Map\'];
%     OutputFolder = [Image_Folder(1:end-1) '-Wave\wave-' num2str(wave_index) '\'];
    map_imgFolder = [Image_Folder(1:end-1) '\RFP_Map\'];
    OutputFolder = [Image_Folder(1:end-1) '\Wave\wave-' num2str(wave_index) '\'];
    WormRegionFolder = [OutputFolder 'worm_region\'];
    WormRegionCheckFolder = [OutputFolder 'worm_region_check\'];

    desired_width = 512;
    desired_height = 512;
    image_format = '.tiff';

    map_imgs = dir([map_imgFolder '*.tiff']);
    b_imgs = dir([WormRegionFolder '*.tiff']);
    image_num = length(b_imgs);
    frame_indexs = zeros(image_num,1);
    for i = 1:image_num
        frame_indexs(i) = str2num(b_imgs(i).name(1:end-5));
    end
    [frame_indexs,sort_seq] = sort(frame_indexs);
    b_imgs = b_imgs(sort_seq);
    Start_Index = frame_indexs(1);
    End_Index = frame_indexs(end);

for i=Start_Index:End_Index
    map_img = imread([map_imgFolder map_imgs(i).name]);
    map_img = imresize(map_img,[desired_height,desired_width]);
    binary_worm_region = imread([WormRegionFolder num2str(i) '.tiff']);
    rgb_img = uint8(zeros(desired_height,desired_width,3));
    rgb_img(:,:,1) = map_img;
%         rgb_img(:,:,1) = 20;
    rgb_img(:,:,2) = map_img;
    rgb_img(:,:,3) = map_img;
    rgb_img(:,:,2) = rgb_img(:,:,2)+binary_worm_region*(75/255);
    
    imwrite(rgb_img, [WormRegionCheckFolder num2str(i) image_format]);
    disp(['Seg Updated: ' num2str(i-Start_Index+1) '/' num2str(image_num)]);
end

end