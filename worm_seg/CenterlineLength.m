function L = CenterlineLength(imgFolder,wave_index)
    frame_rate = 8;
%     waveFolder = [imgFolder 'Wave\'];
    waveFolder = [imgFolder 'Wave\wave-' num2str(wave_index) '\'];
%     wave_peaks = load([imgFolder(1:end-1) '.txt']);
%     wave_peak = wave_peaks(wave_index);
    centerlineFolder = [waveFolder 'centerline\'];
    
    centerlines = dir([centerlineFolder '*.mat']);
    image_num = length(centerlines);
    for i = 1:image_num
        frame_indexs(i) = str2num(centerlines(i).name(1:end-4));
    end
    [~,sort_seq] = sort(frame_indexs);
    centerlines = centerlines(sort_seq);
    
    L = zeros(1,image_num);
    start_index = str2num(centerlines(1).name(1:end-4));
%     load([centerlineFolder centerlines(1).name(1:end-4)]);
    
    for i = 1:image_num
        load([centerlineFolder num2str(start_index+i-1) '.mat']);
        L(i) = sum(sqrt(sum((centerline(1:end-1,:)-centerline(2:end,:)).^2,2)));
    end
       
%     figure;plot((1:image_num)/frame_rate,L);xlabel('Time(s)');ylabel('Worm Length');
%     hold on; line([wave_peak-start_index+1,wave_peak-start_index+1]/frame_rate,[L(wave_peak-start_index+1)-20,L(wave_peak-start_index+1)+20],'color','r');
%     title(['Wave ' num2str(wave_index)])
%     saveas(gcf,[waveFolder 'wormLength.fig']);


end
    