function mvMapImg(Folder, wave_index,channel)
% clear;
% Folder = 'H:\Backup\20201210\F3\';
% wave_index = 1;
% channel = 'r';

% close all
frame_rate = 8;
waveFolder = [Folder 'Wave\wave-' num2str(wave_index) '\'];
load([waveFolder 'frame_seq.mat']);
if strcmp(channel,'r')
    MapFolder = [Folder 'RFP_Map\'];
    WaveMapFolder = [waveFolder 'RFP_Map\'];
    frame_seq = rfp_frame_seq;
elseif strcmp(channel,'g')
    MapFolder = [Folder 'GCaMP_Map\'];
    WaveMapFolder = [waveFolder 'GCaMP_Map\'];
    frame_seq = gcamp_frame_seq;
end

map_imgs = dir([MapFolder,'*.tiff']);

% waveFolder = [Folder(1:end-1) '-Wave\wave-' num2str(wave_index) '\'];
centerlineFolder = [waveFolder 'centerline\'];

if ~exist(WaveMapFolder,'dir')
    mkdir(WaveMapFolder);
end
centerlines = dir([centerlineFolder '*.mat']);
if strcmp(channel,'r')
    for i = frame_seq
        copyfile([MapFolder map_imgs(i).name],[WaveMapFolder num2str(i) '.tiff']);
        disp(['copy: ' num2str(i) '  ' num2str(i-frame_seq(1)+1) '/' num2str(length(frame_seq))]);
    end
elseif strcmp(channel,'g')
    for i = frame_seq
        copyfile([MapFolder map_imgs(i).name],[WaveMapFolder num2str(i) '.tiff']);
        disp(['copy: ' num2str(i) '  ' num2str(i-frame_seq(1)+1) '/' num2str(length(frame_seq))]);
    end
end
end