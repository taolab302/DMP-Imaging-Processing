function neuron_I = SingleNeuronIntensity(neuron_name,wave_index,cali)
    
    Folder = 'G:\Backup\20191129\Intestine-9\';
%     neuron_radius = 4;
    if strcmp(neuron_name,'AVL')
        neuron_radius = 10;
    elseif strcmp(neuron_name,'DVB')
        neuron_radius = 4;
    else
        neuron_radius = 4;
    end
%     wave_index = 12;

    midtime = load([Folder(1:end-1) '.txt']);
    time_frame = midtime(wave_index,2);
%     close all
    frame_seq = (midtime(wave_index,1)-time_frame+1):(midtime(wave_index,1)+time_frame);
    neuron_I = zeros(length(frame_seq),1);

    waveFolder = [Folder(1:end-1) '-Wave\wave-' num2str(wave_index) '\'];
    neuron_pos = load([waveFolder,'neuron\',neuron_name,'.txt']);
%     neuron_pos = neuron_pos;
    imgs = dir([Folder,'*.tiff']);  
    
    for i = 1:length(frame_seq)
        frame_index = frame_seq(i);
        worm_Image = imread([Folder imgs(frame_index).name]);
        
        if strcmp(cali, 'cali')
            load('G:\intestine code\fluo_pattern');
            [worm_Image,Background] = CalibrateImage(worm_Image,inverse_GCaMP_pattern,150);
        end
            
        [neuron_I(i), ~] = ExtractFluoEnergyAndBackground(worm_Image,[neuron_pos(i,1), neuron_pos(i,2)],neuron_radius,1.0);
        disp(['frame: ' num2str(i) '/' num2str(length(frame_seq))])
    end
    
    if strcmp(cali, 'cali')
        save([waveFolder 'neuron/' neuron_name '_cali.mat'],'neuron_I');
    else
        save([waveFolder 'neuron/' neuron_name '.mat'],'neuron_I')
    end

end