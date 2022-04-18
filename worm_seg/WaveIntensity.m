function I = WaveIntensity(Folder,wave_index,frame_seq,channel)
    waveFolder = [Folder 'Wave\wave-' num2str(wave_index) '\'];
    load('H:\intestine code\fluo_pattern');
    if isempty(channel)
        imgFolder = Folder;
        pattern = inverse_GCaMP_pattern;
        waveFolder = [Folder(1:end-1) '-Wave\wave-' num2str(wave_index) '\'];
    elseif strcmp(channel,'r')
        imgFolder = [Folder 'RFP\'];
        pattern = inverse_RFP_pattern;
    elseif strcmp(channel,'g')
        imgFolder = [Folder 'GCaMP\'];
        pattern = inverse_GCaMP_pattern;
    end
        
    image_names = dir([imgFolder '*.tiff']);
    image_num = length(frame_seq);
    I = zeros(image_num,49);

    if ~isempty(channel)
        AVLr_pos = load([waveFolder 'neuron_pos\red\AVL.txt']);
        DVBr_pos = load([waveFolder 'neuron_pos\red\DVB.txt']);
        AVLg_pos = load([waveFolder 'neuron_pos\green\AVL.txt']);
        DVBg_pos = load([waveFolder 'neuron_pos\green\DVB.txt']);
        sync_struc = load([Folder 'sync_struc']);
        sync_struc = sync_struc.sync_struc;
    end
    
    
    for i = 1:image_num

       frame_index = frame_seq(i);        disp(['frame: ' num2str(frame_index)])
       img = imread([imgFolder image_names(frame_index).name]);
       [img,Background] = CalibrateImage(img,pattern,105);
       img = double(imresize(img,[512,512]));
       
       if isempty(channel)
           bw_img = imread([waveFolder 'worm_region\' num2str(frame_index) '.tiff']);
           load([waveFolder 'centerline_inte\' num2str(frame_index) '.mat']);
       elseif strcmp(channel,'r')
           bw_img = imread([waveFolder 'worm_region\' num2str(frame_index) '.tiff']);
           load([waveFolder 'centerline\' num2str(frame_index) '.mat']);
       elseif strcmp(channel,'g')
           load([waveFolder 'centerline_GCaMP_cali\' num2str(frame_index) '.mat']);
           centerline = centerline/4;
           red_img_index = sync_struc.match_index(frame_index);
           red_neuron_index = red_img_index-sync_struc.match_index(frame_seq(1))+1;
           offsetAVL = AVLg_pos(i,[2,1])-AVLr_pos(red_neuron_index,[2,1]);
           offsetDVB = DVBg_pos(i,[2,1])-DVBr_pos(red_neuron_index,[2,1]);
           bw_img = imread([waveFolder 'worm_region\' num2str(red_img_index) '.tiff']);
           [bw_row,bw_col] = find(bw_img>0);
           bw_pt = [bw_row,bw_col]+double(int32((offsetAVL+offsetDVB)/8));       
           [height,width] = size(bw_img);
           bw_img = zeros(height,width);
           bw_img((bw_pt(:,2)-1)*height+bw_pt(:,1)) = 1;  
       end
       
       [left,right,centerline] = ExtractBoundaries(bw_img,centerline);
       for j = 1:49
           x = [left(j,1) left(j+1,1) right(j+1,1) right(j,1)];
           y = [left(j,2) left(j+1,2) right(j+1,2) right(j,2)];
            bw = poly2mask(y,x,512,512);
            I(i,j) = sum(sum(img.*bw))/sum(sum(bw));
       end

            
        
    end



end