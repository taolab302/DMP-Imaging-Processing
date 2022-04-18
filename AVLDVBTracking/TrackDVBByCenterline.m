function TrackDVBByCenterline(Folder, FrameSeq,wave_index,channel) 
% Track DVB by centerline

% wave_index = 1;
waveFolder = [Folder 'Wave\wave-' num2str(wave_index) '\'];
sync_struc = load([Folder 'sync_struc.mat']);
DVB_thres = 300;
DVB_thres_now = DVB_thres;
neuron_area_thres = 3;
search_interval = 8;

if strcmp(channel,'g')
    PosFolder = [waveFolder 'neuron_pos\green\'];
    NeuronFolder = [Folder 'GCaMP_Neuron\'];
    CenterlineFolder = [waveFolder 'centerline_GCaMP\'];
    NeuronFiles = dir([Folder '*.mat']);
elseif strcmp(channel,'r')
    PosFolder = [waveFolder 'neuron_pos\red\'];
    NeuronFolder = [Folder 'RFP_Neuron\'];
    CenterlineFolder = [waveFolder 'centerline_origin\'];
    NeuronFiles = dir([Folder '*.mat']);
end

if ~exist(PosFolder,'dir')
    mkdir(PosFolder);
end

if strcmp(FrameSeq, 'all') == 1
    FrameSeq = 1:length(NeuronFiles);
end
pos = zeros(length(FrameSeq),2);
for ii = 1:length(FrameSeq)
    i = FrameSeq(ii);
    disp(['Processing ' num2str(ii) '/' num2str(length(FrameSeq))]);
    CenterlineFile = load([CenterlineFolder num2str(i) '.mat']);
    
    if strcmp(channel,'r')
        img_name = char(sync_struc.sync_struc.sync_names(find(sync_struc.sync_struc.match_index==i,1),2));
        img = double(imread([Folder 'RFP\' img_name]));
        filter_size = [5,5];
        h = fspecial('gaussian',filter_size,1.5);
        img = imfilter(img,h);
%         img(img<DVB_thres_now) = 0;
        se = strel('disk',10);
        [height,width] = size(img);
        bw_img = imread([waveFolder 'worm_region\' num2str(i) '.tiff']);
        bw_img = double(imresize(bw_img,[height,width]))>0;
        bw_img = imdilate(bw_img,se);
        
        localwidth = 40;
        tail = CenterlineFile.centerline(end,:);
        mask_x = [tail(2)-localwidth, tail(2)+localwidth, tail(2)+localwidth, tail(2)-localwidth];
        mask_y = [tail(1)-localwidth,tail(1)-localwidth, tail(1)+localwidth, tail(1)+localwidth];
        mask_x(mask_x<1) = 1; mask_x(mask_x>width) = width;
        mask_y(mask_y<1) = 1; mask_y(mask_y>height) = height;
        localmask = poly2mask(mask_x, mask_y, height, width);
        localmask_n = localmask.*~bw_img;
        [~,cc_neuron_info] = ExtractAnteriorFrag(localmask_n,img,DVB_thres_now,tail); 
        D_mat = bwdist(~localmask_n);
        neuron_index = find(cc_neuron_info(:,1)>=neuron_area_thres&...
            D_mat((int32(cc_neuron_info(:,4))-1)*height+int32(cc_neuron_info(:,3)))>3&...
            D_mat((int32(cc_neuron_info(:,4))-1)*height+int32(cc_neuron_info(:,3)))<35,1);
        neuron_I = img((int32(cc_neuron_info(:,4))-1)*height+int32(cc_neuron_info(:,3)));
        if ~isempty(neuron_index)
            neuron_index = neuron_index(neuron_I(neuron_index)==max(neuron_I(neuron_index)));
            pos(ii,:) = cc_neuron_info(neuron_index,[3,4]);
            [pos(ii,2),pos(ii,1)] = UpdateNeuronPos(pos(ii,2),pos(ii,1),search_interval,0.3,img);
        else
           figure(64);imagesc(img);axis equal;colormap('gray');caxis([0 300]);title(['Frame ' num2str(i)]);
           n_likely = find(D_mat((int32(cc_neuron_info(:,4))-1)*height+int32(cc_neuron_info(:,3)))>3&...
           D_mat((int32(cc_neuron_info(:,4))-1)*height+int32(cc_neuron_info(:,3)))<25,1);
           
%            n_likely = find(neuron_I==max(neuron_I),1);
           n_likely_area = cc_neuron_info(n_likely,1);
           hold on;plot(cc_neuron_info(:,4),cc_neuron_info(:,3),'r*');
           hold on;plot(cc_neuron_info(n_likely,4),cc_neuron_info(n_likely,3),'g*');hold off;
            load splat
            sound(y,Fs)
            disp(['Cannot find DVB. Area of green point = ' num2str(n_likely_area)]);
            prompt = 'Is the green point right? any number = yes, or choose another pos: ';
            neuron_pos = input(prompt);
            while ~isnumeric(neuron_pos)||(length(neuron_pos)~=2&&length(neuron_pos)~=1)
                prompt = 'Input not numeric! Choose another pos(any number = green point): ';
                neuron_pos = input(prompt);
            end
            if length(neuron_pos)==1
                neuron_pos = cc_neuron_info(n_likely,3:4);
            end
            [neuron_pos(2),neuron_pos(1)] = UpdateNeuronPos(neuron_pos(2),neuron_pos(1),search_interval,0.3,img);
            [DVB_thres_now,~] =  ExtractFluoEnergyAndBackground(img,neuron_pos([2 1]),search_interval,0.3);
            pos(ii,:) = neuron_pos;
        end
        
    elseif strcmp(channel,'g')
        img_name = char(sync_struc.sync_struc.sync_names(i,1));
        NeuronFile = load([NeuronFolder img_name '.mat']);
        neurons = NeuronFile.neurons;
        tail = CenterlineFile.centerline(end-1:end,:);
        dis = inf; best_index = 0;
        for j = 1:size(neurons,1)
            if dist(tail(2,:),neurons(j,:)') < dis && dot((tail(2,:)-tail(1,:)),(neurons(j,:)-tail(2,:))) > 0
                dis = dist(tail(2,:),neurons(j,:)');
                best_index = j;
            end
        end
        if best_index == 0
            pos(ii,:) = tail(2,:);
        else
            pos(ii,:) = neurons(best_index,:);
        end
    end

end
output_name = [PosFolder 'DVB.txt'];
fid = fopen(output_name,'a');  
for i = 1:length(FrameSeq)
    fprintf(fid,'%d    %d\n', pos(i,2), pos(i,1));
end
fclose(fid);

    



        