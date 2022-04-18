function MapGCaMP2RFP_single(Image_Folder, neuron_name,wave_index, rfp_frames,gcamp_frames)
% Map red neuron positions to green images

frame_rate = 16;
% neuron_num = length(tracking_index);

% Tracking Parameters Setting
intensity_ratio = 0.5;
Neuron_Motion_Threshold = 0;
neuron_num = 1;

waveFolder = [Image_Folder 'Wave\wave-' num2str(wave_index) '\'];
posFolder = [waveFolder 'neuron_pos\'];
if ~exist([posFolder 'red\'],'dir')
    mkdir([posFolder 'red\']);
end
GCaMP_Folder = [Image_Folder 'GCaMP\'];
RFP_Folder = [Image_Folder 'RFP\'];
image_num = length(rfp_frames);
sync_res = load([Image_Folder 'sync_struc.mat']);
gcamp_neuron_pos = load([posFolder 'green\' neuron_name '.txt']);

rfp_imgs = dir([RFP_Folder '*.tiff']);

neuron_radius = 8;  % see/record in neuron_radius file

rfp_start = rfp_frames(1);
rfp_end = rfp_frames(image_num);
rfp_neurons_x = zeros(rfp_end-rfp_start+1,neuron_num);
rfp_neurons_y = zeros(rfp_end-rfp_start+1,neuron_num);
Ax_tempt = nan(1,10); Ay_tempt = nan(1,10);

% map neuron positions
for iFrame=rfp_start:rfp_end
    worm_img = imread([RFP_Folder rfp_imgs(iFrame).name]);
    gcamp_index = find(sync_res.sync_struc.match_index==iFrame,1);
    if isempty(find(gcamp_frames==gcamp_index, 1))
        gcamp_index = find(sync_res.sync_struc.match_index==iFrame,1,'last');
    end
    motion_offset = [0,0];
    pos_index = gcamp_index-gcamp_frames(1)+1;

%         [Ax_tempt(1),Ay_tempt(1)] = UpdateNeuronPos(neuron_pos_x,neuron_pos_y,neuron_radius(n),intensity_ratio,worm_img);
        [Ax_tempt(1),Ay_tempt(1)] = UpdateNeuronPos(gcamp_neuron_pos(pos_index,1),gcamp_neuron_pos(pos_index,2),neuron_radius,intensity_ratio,worm_img);
        
        % Update 
        k = 1;
        rfp_neurons_x(iFrame-rfp_start+1,1) = Ax_tempt(k);
        rfp_neurons_y(iFrame-rfp_start+1,1) = Ay_tempt(k);
        disp(['RFP Neuron Pos ' num2str(iFrame) '  Ax = ',num2str(Ax_tempt(k)) ' Ay = ',num2str(Ay_tempt(k))]);
end

% Write all neuron positions into file
if length(rfp_frames) > 1
    output_name = [posFolder 'red\' neuron_name '.txt'];
    fid = fopen(output_name,'w');  
    for i = 1:length(rfp_frames)
        fprintf(fid,'%d    %d\n', rfp_neurons_x(i,1), rfp_neurons_y(i,1));
    end
    fclose(fid);
end

end