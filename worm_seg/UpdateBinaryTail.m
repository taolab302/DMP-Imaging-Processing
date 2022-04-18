% update the tail part of binary worm region, in RFP
Folder = 'H:\Backup\20201210\F2\';
waveFolder = [Folder 'Wave\'];
imgFolder = [Folder 'RFP\'];
map_imgFolder = [imgFolder(1:end-1) '_Map\'];
regionFolder = [waveFolder 'worm_region_test\'];
regioncheckFolder = [waveFolder 'worm_region_check_test\'];
if ~exist(regionFolder,'dir')
    mkdir(regionFolder);
end
if ~exist(regioncheckFolder,'dir')
    mkdir(regioncheckFolder);
end

image_format = '.tiff';
img_names = dir([imgFolder '*.tiff']);

WormRegionPos = load([waveFolder 'WormRegionPos.mat']);
centerline_start = WormRegionPos.centerline_start;
centerline_start = int32(centerline_start);

frame_seq = 1:2879;
Start_Index = frame_seq(1);
End_Index = frame_seq(end);

desired_width = 512; desired_height = 512;
height = 512; width = 512;
localwidth = 20;
 neighbor = [0 0; 0 -1; -1 -1; -1 0; -1 1; 0 1; 1 1; 1 0; 1 -1];

for i = Start_Index:End_Index
    
    binary_worm = double(imread([waveFolder 'worm_region\' num2str(i) '.tiff']));
    img = double(imread([imgFolder img_names(i).name]));
    img = imresize(img,[height,width]);
    
    load([waveFolder 'centerline\' num2str(i) '.mat']);
    tail = centerline(end,:);
    mask_x = [tail(2)-localwidth, tail(2)+localwidth, tail(2)+localwidth, tail(2)-localwidth];
    mask_y = [tail(1)-localwidth, tail(1)-localwidth, tail(1)+localwidth, tail(1)+localwidth];
    mask_x(mask_x<1) = 1; mask_x(mask_x>width) = width;
    mask_y(mask_y<1) = 1; mask_y(mask_y>height) = height;
    localmask = poly2mask(mask_x, mask_y, height, width);
    
    local_img = img.*localmask;
    frag_orig = binary_worm.*localmask;
    
    frag_new = activecontour(local_img, frag_orig);
    
    binary_worm_new = (binary_worm - frag_orig + frag_new)>0;
    binary_worm_new = bwmorph(binary_worm_new,'bridge');
    binary_worm_new = imfill(binary_worm_new,'holes');
    
    map_img = imread([map_imgFolder img_names(i).name]);
    map_img = imresize(map_img,[desired_width,desired_height]);

    % save the binary worm region
        rgb_img = uint8(zeros(desired_height,desired_width,3));
        rgb_img(:,:,1) = map_img;
        rgb_img(:,:,2) = map_img;
        rgb_img(:,:,3) = map_img;
        rgb_img(:,:,2) = rgb_img(:,:,1)+uint8(75*binary_worm_new);
        for j = 1:9
            rgb_img(centerline_start(i-Start_Index+1,1)+neighbor(j,1),centerline_start(i-Start_Index+1,2)++neighbor(j,2),1) = 255;
            rgb_img(centerline_start(i-Start_Index+1,1)+neighbor(j,1),centerline_start(i-Start_Index+1,2)+neighbor(j,2),2) = 0;
            rgb_img(centerline_start(i-Start_Index+1,1)+neighbor(j,1),centerline_start(i-Start_Index+1,2)+neighbor(j,2),3) = 0;    
        end
        
    imwrite(rgb_img, [regioncheckFolder num2str(i) image_format]);
    imwrite(binary_worm_new*255, [regionFolder num2str(i) image_format]);
    
    disp(['Processed image ' num2str(i) '/' num2str(length(frame_seq))])
end

    load handel
    sound(y,Fs)
    
