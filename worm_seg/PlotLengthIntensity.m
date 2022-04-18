function PlotLengthIntensity(waveFolder,wave_index,frame_rate,intestineFolder)
    figure(128);
    pixel_size = 6.5/4; % unit: um
    % wave_index = 9;
    % imgFolder = ['H:\FluoImages\20191128\Intestine-3\'];
    %     waveFolder = [imgFolder(1:end-1) '_wave_' num2str(wave_index) '\'];
    %     wave_peaks = load([imgFolder(1:end-1) '.txt']);
    %     wave_peak = wave_peaks(wave_index);
    centerlineFolder = [waveFolder '\centerline\'];
    centerlines = dir([centerlineFolder '*.mat']);
    load([waveFolder 'waveIntensity.mat']);
    load([waveFolder 'wormLength.mat']);

    % frame_rate = 8;
    image_num = length(L);
    % frame_indexs = zeros(1,image_num);
    %     for i = 1:image_num
    %         frame_indexs(i) = str2num(centerlines(i).name(1:end-4));
    %     end
    % [~,sort_seq] = sort(frame_indexs);
    % centerlines = centerlines(sort_seq);
    % start_index = str2num(centerlines(1).name(1:end-4));

    ah1 = subplot(2,1,1);
    pos1 = get(ah1,'Position');
    % plot((1:image_num)/frame_rate,smooth(L,'moving',4));
    plot((1:image_num)/frame_rate,L*pixel_size);
    hold on;plot((1:image_num)/frame_rate,smooth(L,0.1,'loess'),'k');
    legend('raw','smoothed')
    % line([wave_peak-start_index+1,wave_peak-start_index+1]/frame_rate,[min(L)-10 max(L)+10],'color','r');
    axis([0 image_num/frame_rate min(L*pixel_size)-10 max(L*pixel_size)+10])
    xlabel('Time(s)');ylabel('Worm Length (um)');
    hold on; 
    title(['Wave ' num2str(wave_index)]);
    
    ah2 = subplot(2,1,2);
    pos2 = get(ah2,'Position');
    pos2(3) = pos1(3);
    % hold on;
    % imagesc([0.5,image_num-0.5]/frame_rate,[L(wave_peak-start_index+1)+20,L(wave_peak-start_index+1)+68],I');
    imagesc([0.5,image_num-0.5]/frame_rate,[0 1],I');colormap('jet');
    colorbar;set(ah2,'Position',pos2);
    % line([wave_peak-start_index+1,wave_peak-start_index+1]/frame_rate,[0,1],'color','r');
    yticks([0 1]); yticklabels({'A','P'});
    title('Intensity')
    xlabel('Time(s)');
    % colorbar('southoutside');
    saveas(gca,[intestineFolder 'Intensity-Length wave_' num2str(wave_index) '.fig'])
    fig = figure(128);
    print(fig,[intestineFolder,'Intensity-Length wave_' num2str(wave_index)], '-djpeg','-r300');

    % axis([0 image_num/frame_rate L(wave_peak-start_index+1)-20 L(wave_peak-start_index+1)+70])
end