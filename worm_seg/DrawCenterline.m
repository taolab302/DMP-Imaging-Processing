function DrawCenterline(Folder,frame_seq)
% Draw centerline in current axes
            Image_Folder = [Folder 'worm_region_check\'];
            Centerline_Folder = [Folder 'centerline\'];
            worm_pos = load([Folder 'WormRegionPos.mat']);
            if ~exist([Folder 'fig\'])
                mkdir([Folder 'fig\']);
            end
            worm_region = worm_pos.worm_regions;
%             frame_rate = app.FRAME_RATE;
%             video_name = [app.Current_Folder 'processed_video'];
            
            % Parameters setting
            image_format = '.tiff';
            image_names = dir([Image_Folder, '*' image_format]);
            Start_Index = frame_seq(1);
            End_Index = frame_seq(end);
            
            line_width = 1.2;  % figure paramters
            marker_size = 12;
            
%             writerObj = VideoWriter([video_name '.mp4'],'MPEG-4'); % make video
%             writerObj.FrameRate = frame_rate;
%             open(writerObj);
            
            figure;
            for i=Start_Index:End_Index
                centerline_name = [Centerline_Folder num2str(i) '.mat'];
                % because skip list, some centerlie file does not exist
                if ~exist(centerline_name, 'file')
                    continue;
                end
                data = load(centerline_name);
                centerline = data.centerline;
                
                % load image 
                image_name = [Image_Folder num2str(i) image_format];
                img = imread(image_name);
                worm_img = img;
%                 worm_img = img(worm_region(i-Start_Index+1,1):worm_region(i-Start_Index+1,2),...
%                     worm_region(i-Start_Index+1,3):worm_region(i-Start_Index+1,4),:);
                
                % draw image and plot centerline
                imshow(worm_img);colormap(gray);axis image;hold on;
                plot(centerline(:,2),centerline(:,1),'b-','LineWidth',line_width);
                plot(centerline(1,2),centerline(1,1),'r.','MarkerSize',marker_size,'LineWidth',line_width);
                hold off;
                title(['Image ' num2str(i)]);
                saveas(gcf,[Folder 'fig\' num2str(i) image_format]);
            
%                 % write this frame into video
%                 current_figure = getframe(gcf);
%                 region = current_figure.cdata;
%                 writeVideo(writerObj,region);
            end

end