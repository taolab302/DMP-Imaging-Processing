function GreenCenterlineCali(Folder, wave_index,gcamp_frame_seq)
    waveFolder = [Folder 'Wave\wave-' num2str(wave_index) '\'];
    AVLr_pos = load([waveFolder 'neuron_pos\red\AVL.txt']);
    DVBr_pos = load([waveFolder 'neuron_pos\red\DVB.txt']);
    AVLg_pos = load([waveFolder 'neuron_pos\green\AVL.txt']);
    DVBg_pos = load([waveFolder 'neuron_pos\green\DVB.txt']);
    
    sync_struc = load([Folder 'sync_struc']);
    sync_struc = sync_struc.sync_struc;
    
    CLFolderG = [waveFolder 'centerline_GCaMP\'];
    OutFolder = [waveFolder 'centerline_GCaMP_cali\'];
    if ~exist(OutFolder,'dir')
        mkdir(OutFolder);
    end    
    for i = gcamp_frame_seq
        red_index = sync_struc.match_index(i);
        red_index = red_index-sync_struc.match_index(gcamp_frame_seq(1))+1;
        offsetAVL = AVLg_pos(i-gcamp_frame_seq(1)+1,[2,1])-AVLr_pos(red_index,[2,1]);
        offsetDVB = DVBg_pos(i-gcamp_frame_seq(1)+1,[2,1])-DVBr_pos(red_index,[2,1]);
        centerline = load([CLFolderG num2str(i) '.mat']);
        centerline = centerline.centerline;
        centerline = centerline + (offsetDVB+offsetAVL)/2;
        output_name = [OutFolder num2str(i) '.mat'];
        save(output_name,'centerline');
        disp(['Processed: frame ' num2str(i) '  ' num2str(i-gcamp_frame_seq(1)+1) '/' num2str(length(gcamp_frame_seq))]);
    end
    
   
end