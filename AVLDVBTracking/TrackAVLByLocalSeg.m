function pos = TrackAVLByLocalSeg(Folder, FrameSeq,start_frame,wave_index,channel) 

% find anchor in RFP

localwidth = 240;
locate_thres = 320;
thresA = locate_thres;
AVL_thres = 220;
% AVL_thres = 150;
AnchorAreaA = 140;
AnchorAreaP = 1800;
Bg_thres = 120;
neuron_radius = 8;
% search_interval = 5;
search_interval = 6;
Intensity_Ratio = 0.5;
neuron_area_thres = 1;
neuron_area_prev = neuron_area_thres+20;
AVL_thres_now = AVL_thres;

se = strel('disk',5);

height = 2048;
width = 2048;

GCaMP_Folder = [Folder 'GCaMP\'];
RFP_Folder = [Folder 'RFP\'];
waveFolder = [Folder 'Wave\wave-' num2str(wave_index) '\'];

if strcmp(channel,'g')
    PosFolder = [waveFolder 'neuron_pos\green\'];
elseif strcmp(channel,'r')
    PosFolder = [waveFolder 'neuron_pos\red\'];
    AVL_thres_now = 120;
    AnchorAreaA = 100;
    thresA = 120;
end
if ~exist(PosFolder,'dir')
    mkdir(PosFolder);
end
NeuronFolder = [Folder 'GCaMP_Neuron\'];
NeuronFiles = dir([NeuronFolder '*.mat']);

output_name = [PosFolder 'AVL.txt'];


% anchor_pos = load([PosFolder 'AVL_anchor.txt']);
if strcmp(FrameSeq, 'all') == 1
    FrameSeq = 1:length(NeuronFiles);
end
pos = zeros(length(FrameSeq),2);

greenImgs = dir([GCaMP_Folder '*.tiff']);
redImgs = dir([RFP_Folder '*.tiff']);

sync_struc = load([Folder 'sync_struc.mat']);
WormRegionPos = load([waveFolder  '\WormRegionPos.mat']);
centerline_start = WormRegionPos.centerline_start;

for ii = 1:length(FrameSeq)
    i = FrameSeq(ii);
%     if length(FrameSeq) == 1
%         disp('Replace neuron position by extraction.');
%     else
%         disp(['Processing ' num2str(ii) '/' num2str(length(FrameSeq))]);
%     end
    if strcmp(channel,'g')
%         load('F:\intestine code\fluo_pattern');
%         pattern = inverse_GCaMP_pattern;
        ratio = 12;
        img = double(imread([GCaMP_Folder greenImgs(i).name]));
        cl_start_index = sync_struc.sync_struc.match_index(i);
        cl_start = centerline_start(cl_start_index-sync_struc.sync_struc.match_index(start_frame)+1,:);
%         [img,~] = CalibrateImage(img,pattern,105);
    elseif strcmp(channel,'r')
        ratio = 15;
        img = double(imread([RFP_Folder redImgs(i).name]));
        cl_start_index = i;
        cl_start = centerline_start(i-start_frame+1,:);
    end
%     [height,width] = size(img);
    

    cl_start = cl_start*4;
    
    binary_worm = double(imread([waveFolder 'worm_region\' num2str(cl_start_index) '.tiff']));
    binary_worm = (imresize(binary_worm,[height, width]))>0;
    binary_worm = imdilate(binary_worm,se);
    
    mask_x = [cl_start(2)-localwidth, cl_start(2)+localwidth, cl_start(2)+localwidth, cl_start(2)-localwidth];
    mask_y = [cl_start(1)-localwidth, cl_start(1)-localwidth, cl_start(1)+localwidth, cl_start(1)+localwidth];
    mask_x(mask_x<1) = 1; mask_x(mask_x>width) = width;
    mask_y(mask_y<1) = 1; mask_y(mask_y>height) = height;
    localmask = poly2mask(mask_x, mask_y, height, width);
    
%     local_img = img.*localmask;
    filter_size = [5,5];
    h = fspecial('gaussian',filter_size,1.5);
    local_img = imfilter(img.*localmask,h);
    local_img(local_img<Bg_thres) = 0;
    
    anchorP_bw = (localmask.*binary_worm)>0;
    D_anchorP = bwdist(anchorP_bw);
    
    localmaskA = localmask; localmaskA(binary_worm>0) = 0;
    [ccA,ccA_info] = ExtractAnteriorFrag(localmaskA,local_img,thresA,cl_start);
    dist_anchorP = D_anchorP((int32(ccA_info(:,4))-1)*height+int32(ccA_info(:,3)));
    anchorA = find(ccA_info(:,1)>AnchorAreaA & ccA_info(:,1)<AnchorAreaA*ratio & ...
        localmaskA((int32(ccA_info(:,4))-1)*height+int32(ccA_info(:,3)))>0 & dist_anchorP>32 & dist_anchorP<320);
    
    centerA = [];
    
    while isempty(anchorA)&&length(centerA)~=2
        SaveNeuronPos(output_name,pos,FrameSeq,start_frame);
        figure(64);imagesc(img);axis equal;colormap('gray');caxis([0 500]);title(['Frame ' num2str(i)]);
        hold on; plot(cl_start(2),cl_start(1),'b*');hold on; plot(ccA_info(:,4),ccA_info(:,3),'r*');
        disp(['Cannot find anterior anchor. Current thres = ' num2str(thresA) '  Max area = ' num2str(max(ccA_info(:,1)))]);
        prompt = 'Choose another thres, or choose a center for it: ';
        centerA = input(prompt);
        while ~isnumeric(centerA)
            prompt = 'Input not numeric! Choose anterior anchor center: ';
            centerA = input(prompt);
        end
        
        if length(centerA) == 2
            thresA = ExtractFluoEnergyAndBackground(img,centerA([2 1]),search_interval*3,1);
        elseif length(centerA) == 1
            thresA = centerA;
        end
        [ccA,ccA_info] = ExtractAnteriorFrag(localmaskA,local_img,thresA,cl_start);
        anchorA = find(ccA_info(:,1)>AnchorAreaA & ccA_info(:,1)<AnchorAreaA*ratio & ...
        localmaskA((int32(ccA_info(:,4))-1)*height+int32(ccA_info(:,3)))>0 & D_anchorP((int32(ccA_info(:,4))-1)*height+int32(ccA_info(:,3)))>32 &...
        D_anchorP((int32(ccA_info(:,4))-1)*height+int32(ccA_info(:,3)))<320,1);
    end
 
    anchorA = anchorA(find(ccA_info(anchorA,1)==max(ccA_info(anchorA,1)),1));
%     localmaskP = localmask; localmaskP(ccA.PixelIdxList{anchorA}) = 0;
%     [ccP,ccP_info] = ExtractAnteriorFrag(localmaskP,local_img,locate_thres,cl_start);    
%     anchorP = find(ccP_info(:,1)>AnchorAreaP,1);

    anchorA_bw = zeros(height,width);anchorA_bw(ccA.PixelIdxList{anchorA}) = 1;
    D_anchorA = bwdist(anchorA_bw);
%     anchorP_bw = zeros(height,width);anchorP_bw(ccP.PixelIdxList{anchorP}) = 1;
    
    localmask_n = (localmask-imdilate(anchorA_bw,se)-imdilate(anchorP_bw,se))>0;
    [cc_neuron,cc_neuron_info] = ExtractAnteriorFrag(localmask_n,local_img,AVL_thres_now,cl_start); 
    distA = D_anchorA((int32(cc_neuron_info(:,4))-1)*height+int32(cc_neuron_info(:,3))); 
    distP = D_anchorP((int32(cc_neuron_info(:,4))-1)*height+int32(cc_neuron_info(:,3))); 
%     neuron_index = find(distA<=50&distA>5&cc_neuron_info(:,2)<=130&cc_neuron_info(:,2)>30&cc_neuron_info(:,1)>neuron_area_thres_now);

    direction = dot(cc_neuron_info(:,3:4)-cl_start,cc_neuron_info(:,3:4)-ccA_info(anchorA,3:4),2);
%     neuron_index = find(distA<=80&distA>20 &distP<=140&distP>=30 &...
%         cc_neuron_info(:,1)>neuron_area_thres&cc_neuron_info(:,1)<neuron_area_thres*35 &...
%         cc_neuron_info(:,2)>20 &direction<=0);
    neuron_index = find(distA<=80&distA>20 &distP<=140&distP>=20 &...  %30
        cc_neuron_info(:,1)>=neuron_area_thres&cc_neuron_info(:,1)<=300 &...
        direction<=0);
    
    if ~isempty(neuron_index)
        score = abs((neuron_area_prev-cc_neuron_info(neuron_index,1)));
        n = neuron_index(score == min(score));
        neuron_pos = cc_neuron_info(n(1),3:4);
        neuron_area_prev = cc_neuron_info(n(1),1);
    else
        SaveNeuronPos(output_name,pos,FrameSeq,start_frame);
        figure(64);imagesc(img);axis equal;colormap('gray');caxis([0 200]);title(['Frame ' num2str(i)]);
        hold on; plot(ccA_info(anchorA,4),ccA_info(anchorA,3),'r*');
        hold on; plot(cl_start(2),cl_start(1),'b*');
%         n_likely = find(distA<=50&distA>5&cc_neuron_info(:,2)<=130&cc_neuron_info(:,2)>30);
        n_likely = find(distA<=80&distA>18 &distP<=140&distP>=30 &cc_neuron_info(:,2)>10 &direction<=0);
%         if ii > 1
%             score = (cc_neuron_info(n_likely,1)-I_prev).^2/I_prev+(distA-distA_prev).^2/distA_prev + (distP-distP_prev).^2/distP_prev;
%         else
%            
%         end
        
        n_likely = find(cc_neuron_info(n_likely,1)==max(cc_neuron_info(n_likely,1),1),1);
        n_likely_area = cc_neuron_info(n_likely,1);
        hold on;plot(cc_neuron_info(n_likely,4),cc_neuron_info(n_likely,3),'g*');hold off;
        load splat
        sound(y*0.15,Fs)
        disp(['Cannot find AVL. Area of green point = ' num2str(n_likely_area)]);
        prompt = 'Is the green point right? any number = yes, or choose another pos: ';
        neuron_pos = input(prompt);
        while ~isnumeric(neuron_pos)||(length(neuron_pos)~=2&&length(neuron_pos)~=1)
            prompt = 'Input not numeric! Choose another pos(any number = green point): ';
            neuron_pos = input(prompt);
        end
        if length(neuron_pos)==1
            neuron_pos = cc_neuron_info(n_likely,3:4);
        end
        [AVL_thres_now,~] =  ExtractFluoEnergyAndBackground(local_img,neuron_pos([2 1]),search_interval+2,1);
    end
        [neuron_pos(1),neuron_pos(2)] = UpdateNeuronPos(neuron_pos(1),neuron_pos(2),search_interval,Intensity_Ratio,img);
        pos(ii,:) = neuron_pos;
        disp(['Processing ' num2str(ii) '/' num2str(length(FrameSeq)) '  x = ' num2str(pos(ii,1)) '  y = ' num2str(pos(ii,2)) '  AVL_thres = ' num2str(AVL_thres_now)]);
    
end

SaveNeuronPos(output_name,pos,FrameSeq,start_frame);

end


function SaveNeuronPos(output_name,pos,FrameSeq,start_frame)
        if exist(output_name,'file')
            pos_existed = load(output_name);
        else
            pos_existed = [0 0];
        end
        pos_existed(FrameSeq-start_frame+1,:) = pos(:,[2,1]);
        fid = fopen(output_name,'w');  
        for i = 1:length(pos_existed(:,1))
            fprintf(fid,'%d    %d\n', pos_existed(i,1), pos_existed(i,2));
        end
        fclose(fid);
end
