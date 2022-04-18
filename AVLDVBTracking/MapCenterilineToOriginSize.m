function MapCenterilineToOriginSize(Folder,wave_index,frame_seq)

waveFolder = [Folder 'Wave\wave-' num2str(wave_index) '\'];
CenterlineFolder = [waveFolder 'centerline\'];
NewCenterlineFolder = [waveFolder 'centerline_origin\'];

if ~exist(NewCenterlineFolder, 'dir')
    mkdir(NewCenterlineFolder);
end
CenterlineFiles = dir([CenterlineFolder '*.mat']);
for i = frame_seq
    disp(['Processing ' num2str(i-frame_seq(1)+1) '/' num2str(length(frame_seq))]);
    load([CenterlineFolder num2str(i) '.mat']);
    centerline = 4 * centerline;
    save([NewCenterlineFolder num2str(i) '.mat'], 'centerline');
end