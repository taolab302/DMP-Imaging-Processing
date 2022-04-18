function centerline_start = Segmentation(test,Image_Folder,frame_seq,Worm_Thres,Worm_Area,OutputFolder,fluo_pattern)
    %test 1128-I3 wave_start: 331,940,1351,1936,2462,2962,3411,4060,4529
    partition_num = 50;
    image_names = dir([Image_Folder '*.tiff']);
    
    if strcmp(frame_seq,'all')
        frame_seq = 1:length(image_names);
    end
    image_num = length(frame_seq);
    boundaries = zeros(partition_num,image_num);

    %segmentation
    if strcmp(test,'test')
            img = double(imread([Image_Folder,image_names(frame_seq(1)).name]));
            img = imresize(img,[512,512]);
            map_imgFolder = [Image_Folder(1:end-1) '_Map\'];
            img_fig = imread([map_imgFolder image_names(frame_seq(1)).name]);
            img_fig = imresize(img_fig,[512,512]);
            [height, width] = size(img);
            
             % 补偿光场分布
%             % 法1：中心30%区域减一常值（由图象决定）
%             [img,~] = AmendCenterI(img,height,width,[243,294]);
            
%             % 法2：拟合的高斯pattern矫正，MaxBackground在图像上算
% % %             load('G:\intestine code\fluo_pattern');
%             [height, width] = size(img);
% % %             [~,MaxBackground] = AmendCenterI(img,height,width,[243,294]); %GCaMP
% % %             [img,~] = CalibrateImage(img,imresize(fluo_pattern.inverse_bg_pattern,[512,512]),MaxBackground);
%             [img,MaxBackground] = AmendCenterI(img,height,width,[177.25,286]); %RFP
%             [img,~] = CalibrateImage(img,imresize(fluo_pattern.inverse_RFP_pattern,[512,512]),MaxBackground);

            
            %% 法3：减没有虫时的背景，非KyIs777用
% %               load(fluo_pattern);
% %               rate = median(img,'all')/median(fluo_pattern.GCaMP_background,'all');
% %               GCaMP_background = fluo_pattern.GCaMP_background*rate;
%               GCaMP_background = fluo_pattern.GCaMP_background;
%               img = medfilt2(img,[5,5]) - imresize(GCaMP_background,[height,width]);
%               img = img.*imresize(fluo_pattern.inverse_bg_pattern,[512 512]);
%               img(img<0) = 0;
         %%   
            worm_region = [1, height, 1, width];
            [binary_image,worm_area,~,new_worm_region] = worm_seg_single(img, Worm_Thres, worm_region, Worm_Area);   
%             binary_image = imclose(binary_image,strel('disk',5));
            img_fig = img_fig(new_worm_region(1):new_worm_region(2),new_worm_region(3):new_worm_region(4));
            rgb_img = uint8(zeros(length(binary_image(:,1)),length(binary_image(1,:)),3));
            rgb_img(:,:,1) = img_fig;  %for amend method 1
            rgb_img(:,:,2) = img_fig+uint8(75*binary_image);
            rgb_img(:,:,3) = img_fig;
%             rgb_img(:,:,1) = uint8(img)/1.5+uint8(200*binary_image);
%             rgb_img(:,:,2) = uint8(img)/1.5;
%             rgb_img(:,:,3) = uint8(img)/1.5;
%             imshow(binary_image);colormap(gray);axis image;
            figure(50);imshow(rgb_img);colormap(gray);axis image;   
            disp(['Worm_Thres = ' num2str(Worm_Thres) '  Worm_Area = ' num2str(worm_area)]);
    else
%         worm_seg(Image_Folder,Worm_Thres,Worm_Area,OutputFolder,frame_seq(1),frame_seq(end),fluo_pattern);
        % KyIs777: V3,  else V1
        centerline_start = worm_seg_V3(Image_Folder,Worm_Thres,Worm_Area,OutputFolder,frame_seq(1),frame_seq(end),fluo_pattern);
    end

end