function IntestineCenterline(Folder,frame_seq,head_pos)
% use worm_cv.exe to caculate centerlines of all binary images
            Calculate_Centerline([Folder 'worm_region'], [Folder 'backbone'], frame_seq(1), frame_seq(end));
            
            % load backbone results and generate centerline in original image space
            Centerline_Folder = [Folder 'centerline\'];
            if ~exist(Centerline_Folder,'dir')
                mkdir(Centerline_Folder)
            end
            if ~exist([Folder 'centerline_txt\'],'dir')
                mkdir([Folder 'centerline_txt\'])
            end
            % worm_pos_data = load([app.Current_Folder 'WormRegionPos.mat']);
            Start_Index = frame_seq(1);
            End_Index = frame_seq(end);
            is_reverse = 0;
            
            WormRegionPos = load([Folder 'WormRegionPos.mat']);
            centerline_start = WormRegionPos.centerline_start;
            
            for i=Start_Index:End_Index
                disp(['Processing ceterline backbone file ' num2str(i)]);
                
                
                
                backbone_name = [Folder 'backbone\backbone_' num2str(i) '.bin'];
                % whether backbone file exists
                if ~exist(backbone_name,'file')
                    continue;
                end
                
                % load backbine file
        	backbone = LoadCenterlineResults(backbone_name);
                if ~backbone.length_error
                    centerline = backbone.current_backbone;               
                else
                    centerline = backbone.last_backbone;
                end

                
                % read configure file and get partition number
                config;
                centerline = spline_fitting_partition(centerline, Partition_Num);
               
                
                % add worm region to centerline
                % centerline = centerline + repmat([worm_region(i-Start_Index+1,1) worm_region(i-Start_Index+1,3)],length(centerline),1);
                
                if i == Start_Index && ~isnan(head_pos(1))
                    head_dist = sum((centerline(1,:) - head_pos).^2,2);
        	    tail_dist = sum((centerline(length(centerline),:) - head_pos).^2,2);
        	    if (tail_dist < head_dist)
        	        is_reverse = 1;
        	    end
                end
                
                if is_reverse
                    % reverse centerline
            	    centerline = centerline(end:-1:1,:);
                end
                
                %reset head point to midpoint of anterior edge
                centerline(1,:) = centerline_start(i-Start_Index+1,:);  
                centerline = centerline([1,3:end],:);
%                 centerline = spline_fitting_partition(centerline, Partition_Num);
                
%                 % elongate the centerline till it reaches the boundary of
%                 % head and tail    
                bimg = imread([Folder 'worm_region\' num2str(i) '.tiff']);
                centerline = ElongateCL(centerline,bimg,[],'t');
                
                
                % save centerline
                centerline_name = [Centerline_Folder num2str(i) '.mat'];
                save(centerline_name,'centerline');
                filename = [Folder 'centerline_txt\' num2str(i) '.txt'];
                dump_centerline(centerline,filename);
            end
%             app.CENTERLINE_CALC_FLAG = 1;
end