function DrawWormSeg(Folder,wave_seq)
% Draw centerline in current axes
    for j = 1:length(wave_seq)
        wave_index = wave_seq(j);
        waveFolder = [Folder(1:end-1) '-Wave\wave-' num2str(wave_index) '\'];
            Image_Folder = [waveFolder '\worm_region_check\'];
            Centerline_Folder = [waveFolder '\centerline\'];

%             frame_rate = app.FRAME_RATE;
%             video_name = [app.Current_Folder 'processed_video'];
            
            % Parameters setting
            image_format = '.tiff';
            image_names = dir([Image_Folder, '*' image_format]);
            image_num = length(image_names);
            frame_indexs = zeros(image_num,1);
            for k = 1:image_num
                frame_indexs(k) = str2num(image_names(k).name(1:(end-length(image_format))));
            end
            [frame_indexs,~] = sort(frame_indexs);
            Start_Index = frame_indexs(1);
%             End_Index = frame_seq(end);
%             
%             line_width = 1.2;  % figure paramters
%             marker_size = 12;

            img = imread([Image_Folder num2str(Start_Index) '.tiff']);
            load([Centerline_Folder num2str(Start_Index) '.mat'],'centerline');
            figure(128);clf;imshow(img);axis equal
            hold on;plot(centerline(:,2),centerline(:,1),'bo','markerfacecolor','b','markersize',3)
            for i = 1:49
            dire = centerline(i+1,:) - centerline(i,:);% 从尾部向前搜索
            dire = dire/norm(dire);
            dire_p = [-dire(2) dire(1)];% 中心线法线
            text(centerline(i,2)+dire_p(2)*3,centerline(i,1)+dire_p(1)*3,num2str(i),'color','w','fontsize',8);hold on
            end
            h = figure(128);
            print(h,[Folder(1:end-1) '-Wave\wave' num2str(wave_index) ' seg'], '-djpeg','-r300');
    end
end