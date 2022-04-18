function centerline_start = worm_seg_V3(Image_Folder,Worm_Thres,Worm_Area,OutputFolder,Start_Index,End_Index,fluo_pattern)
    % segment worm region and determine whether skip this image by comparing
    % the area with the first image

    config;
    
    correct_mode = 1;
    channel = 'r';
    
%     local_front_low = 1020;
    local_front_low = 350;
%     local_front_low = 300;
    local_area_low = 150;     
    local_area_high = 1200;
    
%     prompt = 'Choose a worm width: ';
%     worm_width = input(prompt);
    width_thres = 10;
    worm_width = 17;

        neuron_thres = 170;
    if correct_mode ==  1
        neuron_thres = 1000;
    end

    neuron_radius = 10;
    neuron_dist_thres = 13;
    
    if strcmp(channel,'r')
        sync_res = load([Image_Folder(1:end-4) 'sync_struc.mat']);
        map_imgFolder = [Image_Folder(1:end-1) '_Map\'];    % 8-bit images
        green_mapFolder = [Image_Folder(1:end-4) 'GCaMP\'];
    elseif strcmp(channel,'g')
        map_imgFolder = [Image_Folder(1:end-1) '_Map\'];    % 8-bit images
        green_mapFolder = [Image_Folder(1:end-1) '\'];
    end
    
    centerline_start = zeros(length(Start_Index:End_Index),2);

    BackboneFolder = [OutputFolder 'backbone\'];
    WormRegionFolder = [OutputFolder 'worm_region\'];
    WormRegionCheckFolder = [OutputFolder 'worm_region_check\'];
    tempFolder = [OutputFolder 'temp\'];
    if ~exist(WormRegionCheckFolder,'dir')
        mkdir(WormRegionCheckFolder);
    end
    if ~exist(WormRegionFolder,'dir')
        mkdir(WormRegionFolder);
    end
    if ~exist(BackboneFolder,'dir')
        mkdir(BackboneFolder);
    end
    if ~exist(tempFolder,'dir')
        mkdir(tempFolder);
    end

    % desired_width = 2048;
    % desired_height = 2048;
    desired_width = 512;    desired_height = 512;
    step = 2; % for Cali method 1&2
%     tolerence = 0.1;
%     step = 1; % for Cali method 3
    neighbor = [0 0; 0 -1; -1 -1; -1 0; -1 1; 0 1; 1 1; 1 0; 1 -1];

    image_format = '.tiff';
    image_names = dir([Image_Folder, '*' image_format]);
    image_names_green = dir([green_mapFolder, '*' image_format]);

    if Start_Index==0&&End_Index==0
        Start_Index = 0;
        End_Index = length(image_names)-1;
    end

    Skip_List= zeros(length(image_names),1);
    Skip_List_Index = 0;
    
    Worm_Area_prev = Worm_Area;

    if ~exist([OutputFolder 'WormRegionPos.mat'],'file')
        worm_pos = zeros(length(image_names),2);
        worm_regions = zeros(length(image_names),4);
    else
        WormRegionPos = load([OutputFolder 'WormRegionPos.mat']);
        worm_pos = WormRegionPos.worm_pos;
        worm_regions = WormRegionPos.worm_regions;
    end


    for i=Start_Index:End_Index
        if strcmp(channel,'r')
            green_index = find(sync_res.sync_struc.match_index==i,1);
        else
            green_index = i;
        end
        img = double(imread([Image_Folder image_names(i).name]));
        img = imresize(img,[desired_width,desired_height]);
        img_green = double(imread([green_mapFolder image_names_green(green_index).name]));
        img_green = imresize(img_green,[desired_width,desired_height]);

        map_img = imread([map_imgFolder image_names(i).name]);
        map_img = imresize(map_img,[desired_width,desired_height]);
        if i == Start_Index
            [image_height, image_width] = size(img);
            if image_height < desired_height || image_width < desired_width
                disp('Desired height/width is invalid');
                return;
            end
            worm_region = [1,image_height,1,image_width];
        end

        %% Amend Light Ditribution
    %     % Method 1: Deduct a constant at the 30% area of center
    %     [img,~] = AmendCenterI(img,image_height,image_width,[243,294]);

%         % Method 2: Gaussian pattern，MaxBackground calculated in current img
%     %     load('G:\intestine code\fluo_pattern');
% %         [~,MaxBackground] = AmendCenterI(img,image_height,image_width,[243,294]); %GCaMP
% %         [img,~] = CalibrateImage(img,imresize(fluo_pattern.inverse_bg_pattern,[512,512]),MaxBackground);
%         [~,MaxBackground] = AmendCenterI(img,image_height,image_width,[177.25,286]); %RFP
%         [img,~] = CalibrateImage(img,imresize(fluo_pattern.inverse_RFP_pattern,[512,512]),MaxBackground);

        %Method 3: Deduct background in imgs without worm
%         load(fluo_pattern);
%         rate = median(img,'all')/median(fluo_pattern.GCaMP_background,'all');
%         GCaMP_background = fluo_pattern.GCaMP_background*rate;
%         img = medfilt2(img,[5,5]) - imresize(GCaMP_background,[image_height,image_width]);
%         img = img.*imresize(fluo_pattern.inverse_bg_pattern,[512 512]);
%         img(img<0) = 0;

        %% preliminary segment
       [binary_worm_region, Worm_Area, pos, new_worm_region] = worm_seg_single(img, Worm_Thres, worm_region, Worm_Area);
       worm_region_bg = img(new_worm_region(1):new_worm_region(2),new_worm_region(3):new_worm_region(4));
       img_green = img_green(new_worm_region(1):new_worm_region(2),new_worm_region(3):new_worm_region(4));
       
       if i == Start_Index
           Init_Worm_Area = Worm_Area;
           prompt = 'What is the frag center?';
           frag_center_prev = input(prompt);
           prompt = 'CL (head towards) to VD? clockwise = 1; else = -1: ';
           ventral_dir = input(prompt);
       end
       thres_flag = 0;
       if sum(imfill(binary_worm_region,'holes'),'all')>1.3*sum(binary_worm_region,'all')
           thres_flag = Worm_Thres;
           Skip_List = Skip_List(1:Skip_List_Index);
                    raw_worm_pos = worm_pos;
                    load([OutputFolder 'frame_seq.mat']);
                    if exist([OutputFolder 'WormRegionPos.mat'],'file')
                %         worm_pos_info = load([OutputFolder 'WormRegionPos.mat']);                
                        WormRegionPos.centerline_start((Start_Index:End_Index) - rfp_frame_seq(1)+1,:) = centerline_start;
                        centerline_start = WormRegionPos.centerline_start;
                    end
                    % write centerline_start into txt
                    file = fopen([OutputFolder 'centerline_start_xy.txt'],'wt');
                    for k=1:length(centerline_start)
                        fprintf(file,'%d    %d\n', centerline_start(k,2), centerline_start(k,1));
                    end
                    fclose(file);
                    % save worm regions and positions
                        save([OutputFolder 'WormRegionPos.mat'],'worm_pos','worm_regions','raw_worm_pos','Skip_List','centerline_start');
                    % write skip list into backbone folder
                    file = fopen([OutputFolder 'backbone\skiplist.txt'],'wt');
                        for k=1:length(Skip_List)
                            fprintf(file,'%d\n',Skip_List(k));
                        end
                    fclose(file);
                    centerline_start = centerline_start((Start_Index:End_Index) - rfp_frame_seq(1)+1,:);
       end
       while thres_flag >0
           if thres_flag == 0
               worm_region_bg = img(new_worm_region(1):new_worm_region(2),new_worm_region(3):new_worm_region(4));
                rgb_img = uint8(zeros(length(binary_worm_region(:,1)),length(binary_worm_region(1,:)),3));
                rgb_img(:,:,1) = worm_region_bg;  %for amend method 1
                rgb_img(:,:,2) = worm_region_bg+75*binary_worm_region;
                rgb_img(:,:,3) = worm_region_bg;
                figure(50);imshow(rgb_img);colormap(gray);axis image; 
                disp(['Worm_Thres = ' num2str(Worm_Thres) '  Worm_Area = ' num2str(Worm_Area)]);
               prompt = 'Head sticks with tail! Choose a temporary Worm_Thres: ';
               
               Worm_Thres_temp  = input(prompt);
           else
               Worm_Thres_temp = thres_flag;
           end
            [binary_worm_region, Worm_Area, pos,new_worm_region] = worm_seg_single(img, Worm_Thres_temp, worm_region, Worm_Area); 
            worm_region_bg = img(new_worm_region(1):new_worm_region(2),new_worm_region(3):new_worm_region(4));
            rgb_img = uint8(zeros(length(binary_worm_region(:,1)),length(binary_worm_region(1,:)),3));
            rgb_img(:,:,1) = worm_region_bg;  %for amend method 1
            rgb_img(:,:,2) = worm_region_bg+75*binary_worm_region;
            rgb_img(:,:,3) = worm_region_bg;
            figure(50);imshow(rgb_img);colormap(gray);axis image;   
            disp(['Worm_Thres = ' num2str(Worm_Thres_temp) '  Worm_Area = ' num2str(Worm_Area)]);
            prompt = 'Worm_Thres feasible? -1 = yes, or input another value to test: ';
            thres_flag = input(prompt);
       end
       
            
        %% Adapt Worm_Thres using area and the img
            if abs(Init_Worm_Area - Worm_Area) > Init_Worm_Area*Frame_Skip_Thres
                disp(['Worm Lost!  Current Worm_Thres = ' num2str(Worm_Thres) '  Worm_Area = ' num2str(Worm_Area) ' Init_Worm_Area = ' num2str(Init_Worm_Area)]);
                break;
            else
                height = length(binary_worm_region(:,1));
                width = length(binary_worm_region(1,:));
                se = strel('disk',5);
                
                %% Get Anterior Fragment
                 body_cl = bwskel(imerode(binary_worm_region,strel('disk',3)));   % 通过骨架端点选出最接近头的局部区域
                 end_pt_bw = bwmorph(body_cl,'endpoints');
                 end_pt = zeros(sum(end_pt_bw,'all'),2); 
                 [end_pt(:,1),end_pt(:,2)] = find(end_pt_bw>0);
                 dist_frag = sqrt(sum((end_pt-frag_center_prev).^2,2));
                 end_pt = end_pt(dist_frag==min(dist_frag),:);
                 localwidth = 70;
                 mask_x = [end_pt(2)-localwidth, end_pt(2)+localwidth, end_pt(2)+localwidth, end_pt(2)-localwidth];
                 mask_y = [end_pt(1)-localwidth, end_pt(1)-localwidth, end_pt(1)+localwidth, end_pt(1)+localwidth];
                 mask_x(mask_x<1) = 1; mask_x(mask_x>width) = width;
                 mask_y(mask_y<1) = 1; mask_y(mask_y>height) = height;
                 localmask = poly2mask(mask_x, mask_y, height, width);
                 local_I = sum(localmask.*worm_region_bg,'all')/localwidth^2;
                 
%                  if i == Start_Index
%                      local_front_low_prev = local_front_low;
%                      local_I_prev = local_I;
% %                     local_thres_offset = 0;
%                  else
%                     local_front_low = local_front_low_prev + local_I-local_I_prev;             
%                  end
                 
                
                [cc,cc_info] = ExtractAnteriorFrag(binary_worm_region.*localmask,worm_region_bg,local_front_low,frag_center_prev);
             
                cc_info(cc_info(:,1)<local_area_low|cc_info(:,1)>local_area_high,2) = nan;
                temp = localmask.*worm_region_bg;
                for j = 1:length(cc_info(:,1))
                    fragment = cc.PixelIdxList{j};
                    fragment_row = mod(fragment, height);
                    fragment_column = ceil(fragment / height);
                    if range(fragment_row)>50||range(fragment_column)>50||temp(int32(cc_info(j,3)),int32(cc_info(j,4)))==0
                        cc_info(j,2) = nan;
                    end
                end
                frag_index = find(cc_info(:,2)==min(cc_info(:,2)));
                centerA = [];
                
                while isempty(frag_index)&&length(centerA)~=2
                    Skip_List = Skip_List(1:Skip_List_Index);
                    raw_worm_pos = worm_pos;
                    load([OutputFolder 'frame_seq.mat']);
                    if exist([OutputFolder 'WormRegionPos.mat'],'file')
                %         worm_pos_info = load([OutputFolder 'WormRegionPos.mat']);                
                        WormRegionPos.centerline_start((Start_Index:End_Index) - rfp_frame_seq(1)+1,:) = centerline_start;
                        centerline_start = WormRegionPos.centerline_start;
                    end
                    % write centerline_start into txt
                    file = fopen([OutputFolder 'centerline_start_xy.txt'],'wt');
                    for k=1:length(centerline_start)
                        fprintf(file,'%d    %d\n', centerline_start(k,2), centerline_start(k,1));
                    end
                    fclose(file);
                    % save worm regions and positions
                        save([OutputFolder 'WormRegionPos.mat'],'worm_pos','worm_regions','raw_worm_pos','Skip_List','centerline_start');
                    % write skip list into backbone folder
                    file = fopen([OutputFolder 'backbone\skiplist.txt'],'wt');
                        for k=1:length(Skip_List)
                            fprintf(file,'%d\n',Skip_List(k));
                        end
                    fclose(file);
                    centerline_start = centerline_start((Start_Index:End_Index) - rfp_frame_seq(1)+1,:);
                    
                    
                    figure(100);imagesc(worm_region_bg);axis equal;colormap('gray');caxis([0 local_front_low]);title(['Frame ' num2str(i)]);
                    hold on; 
                    plot(cc_info(:,4),cc_info(:,3),'r*');
                    plot(cc_info(~isnan(cc_info(j,2)),4),cc_info(~isnan(cc_info(j,2)),3),'g*');
                    hold off;
                    disp(['Cannot find anterior anchor. Current thres = ' num2str(local_front_low) '  Max area = ' num2str(max(cc_info(:,1)))]);
                    prompt = 'Choose another thres, or choose a center for it: ';
                    centerA = input(prompt);
                    while ~isnumeric(centerA)
                        prompt = 'Input not numeric! Choose anterior anchor center: ';
                        centerA = input(prompt);
                    end

                    if length(centerA) == 2
                        local_front_low = ExtractFluoEnergyAndBackground(img,centerA([2 1]),search_interval*3,1);
                    elseif length(centerA) == 1
                        local_front_low = centerA;
                    end
                    [cc,cc_info] = ExtractAnteriorFrag(binary_worm_region,worm_region_bg,local_front_low,frag_center_prev);
                    cc_info(cc_info(:,1)<local_area_low|cc_info(:,1)>local_area_high,2) = nan;
                    frag_index = find(cc_info(:,2)==min(cc_info(:,2)));
                end
                
                local_front_low_prev = local_front_low;
                local_I_prev = local_I;
                
                fragment = cc.PixelIdxList{frag_index};
                frag_bw = zeros(height,width); frag_bw(fragment) = 1;
                frag_bw = imclose(frag_bw,se);
%                 frag_bw = activecontour(worm_region_bg,frag_bw);
%                 fragment_row = mod(fragment, height);
%                 fragment_column = ceil(fragment / height);
                [fragment_row,fragment_column] = find(frag_bw>0);
                frag_center_prev = [mean(fragment_row), mean(fragment_column)];
%                 if cc_info(frag_index,1)>500
%                     pause;
%                 end

                %% deduct anterior small fragment
                pt = [fragment_row,fragment_column];
                corner_index = convhull(pt);
                corner = pt(corner_index,:);
                frag_mask = poly2mask(pt(corner_index,2),pt(corner_index,1),height, width);
                
                frag_mask_dilate = imdilate(frag_mask,se);
                binary_worm_region_old= binary_worm_region;
                rest = binary_worm_region;
                rest(frag_mask_dilate>0) = 0;
%                 rest(frag_mask>0) = 0;
                cc_rest = bwconncomp(rest);
                for m = 1:length(cc_rest.PixelIdxList)
                      frag_rest = cc_rest.PixelIdxList{m};
                      if length(frag_rest)<local_area_high
                           binary_worm_region(frag_rest) = 0;
                      else
                          rest_new = zeros(height,width);
                          rest_new(frag_rest) = 1;
                      end
                end
                binary_body = binary_worm_region;
%                 binary_worm_region(frag_mask>0) = 1;
%                 D_mat = bwdist(rest_new);
                
                %% locate centerline start 
                D_mat = bwdist(imerode(rest_new,strel('disk',4)));
                edge1 = corner(2:end,:)-corner(1:end-1,:); edge_length1 = sqrt(sum(edge1.^2,2));
                edge2 = [corner(end-1,:);corner(1:end-2,:)]-corner(1:end-1,:);edge_length2 = sqrt(sum(edge2.^2,2));
                angles = sum(edge1.*edge2,2)./(edge_length1.*edge_length2);
                corner = corner(angles>-0.95&angles<sqrt(2)/2,:); 
                corner = [corner; corner(1,:)];
                corner_dist = zeros(length(corner(:,1))-1,1);
                edge_length = zeros(length(corner(:,1))-1,1);
                for m = 1:length(corner(:,1))-1
                  midpoint = [mean(corner(m:m+1,1)),mean(corner(m:m+1,2))];
            %       direction = midpoint - frag_center;
            %       direction = direction/norm(direction);
                  edge_length(m) = norm(corner(m,:)-corner(m+1,:));
                  if norm(corner(m,:)-corner(m+1,:))>worm_width/3
                      corner_dist(m) = D_mat(int32(midpoint(1)),int32(midpoint(2)));
                  else
                      corner_dist(m) = -D_mat(int32(midpoint(1)),int32(midpoint(2)));
                  end
            %      corner_dist(m) = min(pt(corner_index(m),:)-) 
                end
               sec_max = max(corner_dist(corner_dist~=max(corner_dist)));
               if abs(sec_max-max(corner_dist))<=1
                   m1 = find(corner_dist == max(corner_dist)&corner_dist>0);
                   m2 = find(corner_dist == sec_max&corner_dist>0);
                   if edge_length(m1)>edge_length(m2)
                       m = m1;
                   else
                       m = m2;
                   end
               else
                   m = find(corner_dist == max(corner_dist)&corner_dist>0);
               end
               midpoint = [mean(corner(m:m+1,1)),mean(corner(m:m+1,2))];
               
             %% split body edge and identify head & tail
                body_edge = edge(binary_body);
                body_edge = bwmorph(body_edge,'bridge');
                
                edge_img = imfill(body_edge,'holes');
                body_cl = bwskel(imerode(binary_body,strel('disk',3)));
                branch_bw = bwmorph(body_cl,'branchpoints');
                [branch_pty,branch_ptx] = find(branch_bw>0);
                branch_pt = [branch_pty,branch_ptx];
                branch_local = zeros(length(branch_ptx)*9,2);
                for k = 1:length(branch_ptx)
                    branch_local((k-1)*9+1:k*9,:) = branch_pt(k,:)+neighbor;
                end
                body_cl((branch_local(:,2)-1)*height+branch_local(:,1)) = 0;
                body_cl = bwmorph((bwareaopen(body_cl,10)+branch_bw)>0,'bridge');

                end_pt_bw = bwmorph(body_cl,'endpoints');
                end_pt = zeros(2,2); 
                if i==791
                    test = 0;
                end
                [end_pt(:,1),end_pt(:,2)] = find(end_pt_bw>0);
                end_dist = sqrt(sum((end_pt - frag_center_prev).^2,2));
                head_end = end_pt(find(end_dist == min(end_dist),1),:); 
                body_cl_pt = GetCurvePt(body_cl,head_end);
                body_cl_pt = spline_fitting_partition(body_cl_pt, Partition_Num);

                head_dir = midpoint - body_cl_pt(2,:);head_dir = head_dir/norm(head_dir); 
                tail_dir = body_cl_pt(end,:) - body_cl_pt(end-1,:);tail_dir = tail_dir/norm(tail_dir);  
                for j=0:width        
                    head_edge = int32(j*head_dir + body_cl_pt(1,:));                
                    if edge_img(head_edge(1),head_edge(2))==0
                        head_edge = int32((j-1)*head_dir + body_cl_pt(1,:));              
                        break;
                    end       
                end
                for j=0:width        
                    tail_edge = int32(j*tail_dir + body_cl_pt(end,:));                
                    if edge_img(tail_edge(1),tail_edge(2))==0
                        tail_edge = int32((j-1)*tail_dir + body_cl_pt(end,:));
                        break;
                    end       
                end           
                edge_pt = GetCurvePt(body_edge,double(head_edge));
                tail_edge_dist = sqrt(sum((double(tail_edge)-edge_pt).^2,2));
                tail_edge_index = find(tail_edge_dist==min(tail_edge_dist),1);
%                 edge1 = zeros(tail_edge_index,2); edge1_dist = zeros(tail_edge_index,1);
%                 edge2 = zeros(length(edge_pt(:,1))-tail_edge_index,2);edge2_dist = zeros(length(edge_pt(:,1))-tail_edge_index,1);
                edge1 = edge_pt(1:tail_edge_index,:); 
                edge2 = edge_pt(tail_edge_index:end,:);
         
                %% identify ventral/dorsal (ventral/dorsal) edge and make up the ventral edge
                
                edge1_mid = edge_pt(ceil(tail_edge_index/10),:);
                edge2_mid = edge_pt(floor(length(edge_pt(:,1))-tail_edge_index/10),:);
                edge12 = edge2_mid-edge1_mid;
                edge12_mid = (edge1_mid+edge2_mid)/2;
                center_dir = double(head_edge) - edge12_mid;             
                if (center_dir(1)*edge12(2)-center_dir(2)*edge12(1))*ventral_dir>0
                    ventral_edge = edge1;dorsal_edge = edge2;
    %                 ventral_edge_dist = edge1_dist; dorsal_edge_dist = edge2_dist;
                else
                    ventral_edge = edge2(end:(-1):1,:);dorsal_edge = edge1(end:(-1):1,:); %make head the start point
    %                 ventral_edge_dist = edge2_dist(end:(-1):1,:); dorsal_edge_dist = edge1_dist(end:(-1):1,:);
                end

                dorsal_edge_sparse = spline_fitting_partition(dorsal_edge, floor(Partition_Num*0.4));
                dorsal_edge_sparse = dorsal_edge_sparse(2:end-1,:);
                dorsal_edge_sparse = spline_fitting_partition(dorsal_edge_sparse, length(dorsal_edge(:,1))+20);
                dorsal_edge_sparse = unique(int32(dorsal_edge_sparse),'rows','stable');
                dorsal_sparse_bw = zeros(height,width);dorsal_sparse_bw((dorsal_edge_sparse(:,2)-1)*height+dorsal_edge_sparse(:,1)) = 1;
                dorsal_edge_bw = zeros(height,width); dorsal_edge_bw((dorsal_edge(:,2)-1)*height+dorsal_edge(:,1)) = 1;
                ventral_edge_bw = zeros(height,width); ventral_edge_bw((ventral_edge(:,2)-1)*height+ventral_edge(:,1)) = 1;
                [D_dorsal_edge,id_dorsal] = bwdist(dorsal_sparse_bw);
                D_ventral_edge = bwdist(ventral_edge_bw);
%                 dist_dorsal = D_dorsal_edge((ventral_edge(:,2)-1)*height+ventral_edge(:,1));
    %             dist_dorsal_neuron = sum(sum(D_dorsal_edge.*binary_neuron))/sum(sum(binary_neuron));
    
                ventral_edge_length = length(ventral_edge(:,1));
                skip_num = floor(ventral_edge_length*0.21);
%                 skip_num = floor(ventral_edge_length*0.15);
                
                dist_VD = D_dorsal_edge((ventral_edge(:,2)-1)*height+ventral_edge(:,1));
%                 dist_VD = smooth(dist_VD,'moving',4);
                id_dorsal = id_dorsal((ventral_edge(:,2)-1)*height+ventral_edge(:,1));
%                 width_seg = ((1:ventral_edge_length)-skip_num)*(dist_VD(ventral_edge_length-skip_num)-dist_VD(skip_num))/(ventral_edge_length-2*skip_num)+dist_VD(skip_num);
                fitpts_width = floor(skip_num/3);
                fitpts = [(skip_num-fitpts_width:skip_num+fitpts_width)',dist_VD(skip_num-fitpts_width:skip_num+fitpts_width)];
                fitpts = [fitpts;(ventral_edge_length-skip_num-fitpts_width:ventral_edge_length-skip_num+fitpts_width)',dist_VD(ventral_edge_length-skip_num-fitpts_width:ventral_edge_length-skip_num+fitpts_width)];
                p = double(polyfit(fitpts(:,1),fitpts(:,2),1));
                width_seg = (1:ventral_edge_length)*p(1)+p(2);

                new_ventral_edge = ventral_edge;
                for j = ceil(1.2*skip_num):ventral_edge_length-skip_num
                    if dist_VD(j)<width_seg(j)
                       temp = double([mod(id_dorsal(j),height),floor(id_dorsal(j)/height)]);
                       new_ventral_edge(j,:) = (ventral_edge(j,:)-temp)*width_seg(j)/dist_VD(j)+temp;
                    end
                end
                new_ventral_edge = spline_fitting_partition(new_ventral_edge, length(ventral_edge(:,1))*2);
                new_ventral_edge = unique(int32(new_ventral_edge),'rows','stable');
                temp = zeros(height,width);
                temp((new_ventral_edge(:,2)-1)*height+new_ventral_edge(:,1))=1;
                new_ventral_edge = temp;
                
                
                
                
                
%                 mid_body_width = 0.5*(D_dorsal_edge(ventral_edge(skip_num,1),ventral_edge(skip_num,2))+D_dorsal_edge(ventral_edge(ventral_edge_length-2*skip_num,1),ventral_edge(ventral_edge_length-2*skip_num,2)));
%                 new_ventral_edge = MakeUpventralEdge(mid_body_width,ventral_edge,D_dorsal_edge,D_ventral_edge,ventral_edge_bw,skip_num);
%                 binary_worm_region_temp = bwmorph((new_ventral_edge+binary_body+dorsal_edge_bw+ventral_edge_bw+dorsal_sparse_bw)>0,'bridge');
                binary_worm_region_temp = (new_ventral_edge+binary_body+dorsal_edge_bw+ventral_edge_bw+dorsal_sparse_bw)>0;
                binary_worm_region_temp = imfill(binary_worm_region_temp,'holes');
                binary_worm_region_temp = bwareaopen(binary_worm_region_temp,50);
                body_width = CalculateBodyWidth(binary_worm_region_temp,dorsal_edge);
                pass_signal = Inf;
                
%                 %% Quality control
%                 if min(body_width(skip_num:(end-floor(2.5*skip_num))))<width_thres
%                     current_edge = (D_dorsal_edge>=body_width-2&D_dorsal_edge<body_width).*(D_ventral_edge<13);
%                     PlotSegCheck(current_edge,midpoint,i,ventral_edge,dorsal_edge,D_dorsal_edge,binary_worm_region_temp,worm_region_bg);
%                     while pass_signal>0&&pass_signal~=1
%                         if pass_signal == Inf
%                             load splat
%                             sound(y,Fs)
%                         end
%                         if pass_signal>0&&pass_signal~=Inf&&pass_signal~=1
%                             mid_body_width = pass_signal;
%                         else
%                             prompt = ['Current width = ' num2str(mid_body_width) '  Half Width = ' num2str(min(body_width(skip_num:end-skip_num))) '  Choose another width(1=force pass): '];
%                             mid_body_width = input(prompt);
%                             while isempty(mid_body_width)||~isnumeric(mid_body_width)
%                                 prompt = 'Input not a number! Choose a mid-body width again(1=force pass)';
%                                 mid_body_width = input(prompt);
%                             end
%                         end
%                         if mid_body_width~=1
%                             ventral_edge_bw = MakeUpventralEdge(mid_body_width,ventral_edge,D_dorsal_edge,D_ventral_edge,ventral_edge_bw,skip_num);
%                             binary_worm_region_temp = bwmorph((new_ventral_edge+binary_body+dorsal_edge_bw+ventral_edge_bw+dorsal_sparse_bw)>0,'bridge');
%                             binary_worm_region_temp = imfill(binary_worm_region_temp,'holes');
%                             body_width = CalculateBodyWidth(binary_worm_region_temp,dorsal_edge);
%                             current_edge = (D_dorsal_edge>=mid_body_width-2&D_dorsal_edge<mid_body_width).*(D_ventral_edge<13);
%                             PlotSegCheck(current_edge,midpoint,i,ventral_edge,dorsal_edge,D_dorsal_edge,binary_worm_region_temp,worm_region_bg); 
%                             prompt = 'The mid-body width is good?  yes = 0, force pass = 1, or choose another value as width: ';
%                             pass_signal = input(prompt);
%                             while isempty(mid_body_width)||~isnumeric(pass_signal)
%                                 disp('Input not a number! Please retype: ');
%                                 prompt = 'The mid-body width is good?  yes = 0, force pass = 1, or choose another value as width: ';
%                                 pass_signal = input(prompt);
%                             end
%                         else
%                             pass_signal = 1;
%                         end
%                     end
%                 end

                binary_worm_region = binary_worm_region_temp;
%                 
              %% shrink tail to exclude possible DVB
                localwidth = 6;
                mask_x = double([tail_edge(2)-localwidth, tail_edge(2)+localwidth, tail_edge(2)+localwidth, tail_edge(2)-localwidth]);
                mask_y = double([tail_edge(1)-localwidth, tail_edge(1)-localwidth, tail_edge(1)+localwidth, tail_edge(1)+localwidth]);
                mask_x(mask_x<1) = 1; mask_x(mask_x>width) = width;
                mask_y(mask_y<1) = 1; mask_y(mask_y>height) = height;
                localmask = poly2mask(mask_x, mask_y, height, width);
                local_img = worm_region_bg.*localmask;
                frag_orig = binary_worm_region.*localmask;

                frag_new = activecontour(local_img, frag_orig,20);
             
                if sum(frag_new,'all')<sum(frag_orig,'all')
                    binary_worm_region = (binary_worm_region - frag_orig + frag_new)>0;
                    binary_worm_region = bwmorph(binary_worm_region,'bridge');
                    binary_worm_region = imfill(binary_worm_region,'holes');
                end
                
            end
            
            disp(['Processed img ' num2str(i) '  Worm_Thres = ' num2str(Worm_Thres) '  Neuron_thres = ' num2str(neuron_thres) '  Width = ' num2str(median(body_width(skip_num:(end-floor(2.5*skip_num))))) '  Area = ' num2str(Worm_Area)  ' (' num2str(i-Start_Index+1) '/' num2str(End_Index-Start_Index+1) ')']);
            desired_worm_region = GetDesireImageRegion(pos, [desired_height,desired_width],[image_height,image_width]);

        % crop a desired image size (heigth, width)
        if i > Start_Index && abs(Init_Worm_Area - Worm_Area) > Init_Worm_Area*Frame_Skip_Thres
            % the first image must be correct!
            disp(['Skip image: ' num2str(i)]);
            Skip_List_Index = Skip_List_Index + 1;
            Skip_List(Skip_List_Index) = i;
            Worm_Area = Init_Worm_Area;
            worm_pos(i+1,:) = worm_pos(i,:);
        else
            worm_pos(i+1,:) = pos;
        end

        % update worm region and get desired binary image
        worm_region = new_worm_region;
        whole_binary_img = false(size(img));
        whole_binary_img(worm_region(1):worm_region(2),worm_region(3):worm_region(4)) = binary_worm_region;
        binary_worm_region = ...
            whole_binary_img(desired_worm_region(1):desired_worm_region(2),desired_worm_region(3):desired_worm_region(4));
        worm_regions(i+1,:) = desired_worm_region;
%         centerline_start(i-Start_Index+1,:) = [int32(head_edge(1)+worm_region(1)-1),int32(head_edge(2)+worm_region(3)-1)];
        centerline_start(i-Start_Index+1,:) = [int32(midpoint(1)+worm_region(1)-1),int32(midpoint(2)+worm_region(3)-1)];

% %         % deduct long-thin branches
%         binary_worm_region = modify_binary_worm_region(~binary_worm_region,Worm_Area);

        % save the binary worm region
            rgb_img = uint8(zeros(desired_height,desired_width,3));
            rgb_img(:,:,1) = map_img(desired_worm_region(1):desired_worm_region(2),desired_worm_region(3):desired_worm_region(4));
    %         rgb_img(:,:,1) = 20;
            rgb_img(:,:,2) = map_img(desired_worm_region(1):desired_worm_region(2),desired_worm_region(3):desired_worm_region(4));
            rgb_img(:,:,3) = map_img(desired_worm_region(1):desired_worm_region(2),desired_worm_region(3):desired_worm_region(4));
            rgb_img(:,:,2) = rgb_img(:,:,1)+uint8(75*binary_worm_region);
            for j = 1:9
                rgb_img(centerline_start(i-Start_Index+1,1)+neighbor(j,1),centerline_start(i-Start_Index+1,2)++neighbor(j,2),1) = 255;
                rgb_img(centerline_start(i-Start_Index+1,1)+neighbor(j,1),centerline_start(i-Start_Index+1,2)+neighbor(j,2),2) = 0;
                rgb_img(centerline_start(i-Start_Index+1,1)+neighbor(j,1),centerline_start(i-Start_Index+1,2)+neighbor(j,2),3) = 0;    
            end
            
            
            
        imwrite(rgb_img, [WormRegionCheckFolder num2str(i) image_format]);
        imwrite(binary_worm_region*255, [WormRegionFolder num2str(i) image_format]);
    end

    Skip_List = Skip_List(1:Skip_List_Index);
    raw_worm_pos = worm_pos;
    % worm_pos = WormPos_Filtering(worm_pos);
    
    if exist([OutputFolder 'WormRegionPos.mat'],'file')
%         worm_pos_info = load([OutputFolder 'WormRegionPos.mat']);
        load([OutputFolder 'frame_seq.mat']);
        WormRegionPos.centerline_start((Start_Index:End_Index) - rfp_frame_seq(1)+1,:) = centerline_start;
        centerline_start = WormRegionPos.centerline_start;
    end
    
    % write centerline_start into txt
    file = fopen([OutputFolder 'centerline_start_xy.txt'],'wt');
    for i=1:length(centerline_start)
        fprintf(file,'%d    %d\n', centerline_start(i,2), centerline_start(i,1));
    end
    fclose(file);

    % save worm regions and positions
        save([OutputFolder 'WormRegionPos.mat'],'worm_pos','worm_regions','raw_worm_pos','Skip_List','centerline_start');

    
    % write skip list into backbone folder
    file = fopen([OutputFolder 'backbone\skiplist.txt'],'wt');
        for i=1:length(Skip_List)
            fprintf(file,'%d\n',Skip_List(i));
        end
    fclose(file);
    load splat
    sound(y*0.3,Fs)
end

function worm_region = GetDesireImageRegion(worm_pos, desired_size, original_size)
    %% crop new worm region with width and height

    height = desired_size(1);
    width = desired_size(2);
    img_height = original_size(1);
    img_width = original_size(2);

    pos = round(worm_pos);
    row_min = max(1,pos(1)-height/2);
    row_max = min(pos(1) + height/2 - 1, img_height);
    col_min = max(1, pos(2)-width/2);
    col_max = min(img_width, pos(2)+width/2 - 1);

    if (row_max - row_min) ~= height
        if row_max == img_height
            row_min = img_height - height + 1;
        else
            row_max = row_min + height - 1;
        end
    end

    if (col_max - col_min) ~= width
        if col_max == img_width
            col_min = img_width - width + 1;
        else
            col_max = col_min + width - 1;
        end
    end
    worm_region = [row_min, row_max, col_min, col_max];

end

function new_ventral_edge = MakeUpventralEdge(mid_body_width,ventral_edge,D_dorsal_edge,D_ventral_edge,ventral_edge_bw,skip_num)
    [height,width] = size(D_dorsal_edge);
%     new_ventral_edge = zeros(size(ventral_edge_bw));
    
    standby_bw = (D_dorsal_edge>=mid_body_width-2&D_dorsal_edge<=mid_body_width).*(D_ventral_edge<13)>0;
    standby_bw = bwskel(standby_bw);standby_bw([1,height],:) = 0; standby_bw(:,[1,width]) = 0;
    standby_bw = bwareaopen(standby_bw,20);
    standby_endpt = bwmorph(standby_bw,'endpoints');
    [endpt_row,endpt_col] =  find(standby_endpt>0);
    standby_pt = GetCurvePt(standby_bw,[endpt_row(1),endpt_col(1)]);
    
    D_standby = bwdist(standby_bw);
    dist_standby = D_standby((ventral_edge(:,2)-1)*height + ventral_edge(:,1));
    replace_start = 0;replace_end = 0;
    for j = floor(1.2*skip_num):length(ventral_edge(:,1))
        if dist_standby(j)<2
            replace_start = j;break;
        end
    end
    if replace_start == 0
        replace_start = floor(length(ventral_edge(:,1))/3);
    end
    for j = (length(ventral_edge(:,1))-skip_num):(-1):2*skip_num
        if dist_standby(j)<2
            replace_end = j;break;
        end
    end
    if replace_end == 0
        replace_end = ceil(length(ventral_edge(:,1))/3*2);
    end

    new_ventral_edge = ventral_edge_bw;
    exclude_pt = ventral_edge(replace_start:replace_end,:);
    new_ventral_edge((exclude_pt(:,2)-1)*height+exclude_pt(:,1)) = 0; 
    
    dist_replace_start = sqrt(sum((ventral_edge(replace_start,:)-standby_pt).^2,2));
    dist_replace_end = sqrt(sum((ventral_edge(replace_end,:)-standby_pt).^2,2));
    standby_end1 = find(dist_replace_start ==min(dist_replace_start), 1 );
    standby_end2 = find(dist_replace_end ==min(dist_replace_end), 1 );
    if standby_end1<standby_end2
        standby_start = standby_end1;standby_end = standby_end2;
    else
        standby_start = standby_end2;standby_end = standby_end1;
    end
    standby_pt = standby_pt(standby_start:standby_end,:);
    new_ventral_edge((standby_pt(:,2)-1)*height+standby_pt(:,1))=1;
    
end

function body_width = CalculateBodyWidth(binary_worm_region,dorsal_edge)
    height = length(binary_worm_region(:,1));
%     test = bwdist(~binary_worm_region);
    new_cl = bwskel(binary_worm_region);
    new_cl_endpt = bwmorph(new_cl,'endpoints');
    [endpt_row,endpt_col] = find(new_cl_endpt);
%     new_cl_pt = GetCurvePt(new_cl,[endpt_row(1),endpt_col(1)]);
    D_cl = bwdist(new_cl);
    body_width = 2*D_cl((dorsal_edge(:,2)-1)*height+dorsal_edge(:,1));    
end

function PlotSegCheck(current_edge,midpoint,image_index,ventral_edge,dorsal_edge,D_dorsal_edge,binary_worm_region,worm_region_bg)
    figure(200);subplot(1,2,1); 
    current_edge = bwskel(current_edge>0);
    imagesc(D_dorsal_edge+10*current_edge);
    hold on;plot(midpoint(2),midpoint(1),'m*');
    colorbar;axis equal; caxis([0 35]);title(['Image ' num2str(image_index)]);
    hold on;plot(ventral_edge(:,2),ventral_edge(:,1),'r');
    hold on;plot(dorsal_edge(:,2),dorsal_edge(:,1),'g');hold off;
    subplot(1,2,2);imagesc(5*binary_worm_region+10*current_edge+worm_region_bg./150);axis equal;
    hold on;plot(ventral_edge(:,2),ventral_edge(:,1),'r','linewidth',1.5);
    hold on;plot(dorsal_edge(:,2),dorsal_edge(:,1),'g','linewidth',1.5);hold off    
end
