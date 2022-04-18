function neuron_I = SingleNeuronIntensity_RG(neuron_name,channel,cali,Folder, wave_index)
    
%     Folder = 'H:\Backup\20201210\F4\';
    load('..\fluo_pattern\fluo_pattern.mat');
%     wave_index = 4;
    waveFolder = [Folder 'Wave\wave-' num2str(wave_index) '\'];
    load([waveFolder 'frame_seq.mat']);
    
    if strcmp(channel,'r')
        neuron_pos = load([waveFolder,'neuron_pos\red\',neuron_name,'.txt']);
        imgFolder = [Folder 'RFP\'];
        pattern = inverse_RFP_pattern;
        frame_seq = rfp_frame_seq;
    elseif strcmp(channel,'g')
        neuron_pos = load([waveFolder,'neuron_pos\green\',neuron_name,'.txt']);
        imgFolder = [Folder 'GCaMP\'];
        pattern = inverse_GCaMP_pattern;
        frame_seq = gcamp_frame_seq;
    end
    
%     neuron_radius = 4;
    if strcmp(neuron_name,'AVL')
        if strcmp(channel,'g')
            neuron_radius = 10;
        elseif strcmp(channel, 'r')
            neuron_radius = 10;  %9
        end
    elseif strcmp(neuron_name,'DVB')
        if strcmp(channel,'g')
            neuron_radius = 6;  %6
        elseif strcmp(channel, 'r')
            neuron_radius = 6;
        end
    else
        neuron_radius = 4;
    end
%     wave_index = 12;
    imgs = dir([imgFolder,'*.tiff']);
    
    
%     midtime = load([Folder(1:end-1) '.txt']);
%     time_frame = midtime(wave_index,2);
% %     close all
%     frame_seq = (midtime(wave_index,1)-time_frame+1):(midtime(wave_index,1)+time_frame);
    neuron_I = zeros(length(frame_seq),1);

%     waveFolder = [Folder 'Wave\'];

%     neuron_pos = neuron_pos;
      
    
    for i = 1:length(frame_seq)
        frame_index = frame_seq(i);
        worm_Image = imread([imgFolder imgs(frame_index).name]);
        
        if strcmp(cali, 'cali')
            
            [worm_Image,Background] = CalibrateImage(worm_Image,pattern,105);
        end
            
        [neuron_I(i), ~] = ExtractFluoEnergyAndBackground(worm_Image,[neuron_pos(i,1), neuron_pos(i,2)],neuron_radius,1.0);
        disp(['frame: ' num2str(i) '/' num2str(length(frame_seq))])
    end
    
    if strcmp(cali, 'cali')
        save([waveFolder 'neuron_pos\' neuron_name  channel '_cali.mat'],'neuron_I');
    else
        save([waveFolder 'neuron_pos\' neuron_name  channel '.mat'],'neuron_I')
    end

end