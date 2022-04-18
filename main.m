clear;
Folder = 'H:\Backup\20201211\F1\';
wave_index = 1;
rfp_frame_seq = 1:800;
gcamp_frame_seq = 1:1200;  % RFP 16Hz, GCaMP 24Hz

% Segmentation Parameters
Worm_Thres = 120;  % for Cali method 1 & 2
% Worm_Thres = 5;
Worm_Area = 5000;
frame_rate = 16;

pixel_size = 6.5/4; % um
frame_rateR = 16;frame_rateG = 24;
fluo_pattern = load('fluo_pattern\fluo_pattern.mat');

name_elements = strsplit(Folder,'\');
waveFolder = [Folder 'Wave\wave-' num2str(wave_index) '\'];
if ~exist(waveFolder,'dir')
    mkdir(waveFolder);
end
Image_Folder = [Folder 'RFP\'];
%load([waveFolder 'frame_seq.mat']);
frame_seq = rfp_frame_seq;

addpath 'AVLDVBTracking\';
addpath 'worm_seg\';
addpath 'neuron_seg\';
addpath 'common';


%% Worm segmentation (performed in wCherry images)
prompt = 'Execute segmentation? 0 = yes, other values = skip: ';
seg_flag = input(prompt);
if seg_flag ==0
    OutputFolder =  waveFolder;
    if ~exist(OutputFolder,'dir')
        mkdir(OutputFolder);
    end
    Segmentation('test',Image_Folder,frame_seq,Worm_Thres,Worm_Area,OutputFolder,fluo_pattern);
    prompt = 'Worm_Thres feasible? 0 = yes, or input another value to test: ';
    thres_flag = input(prompt);
    while thres_flag>0
        Segmentation('test',Image_Folder,frame_seq,thres_flag,Worm_Area,OutputFolder,fluo_pattern);
        prompt = 'Worm_Thres feasible? 0 = yes, or input another value to test: ';
        Worm_Thres = thres_flag;
        thres_flag = input(prompt);
    end
    tic
    Segmentation(1,Image_Folder,frame_seq,Worm_Thres,Worm_Area,OutputFolder,fluo_pattern);
%     Segmentation(1,Folder,frame_seq(18:end),Worm_Thres,Worm_Area,wave_index);
    toc
end


%% Head position and calculate worm length 
CorrectCLstart(OutputFolder);
prompt = 'What is the head position? ';
head_pos = input(prompt);
% MIJ.closeAllWindows;
IntestineCenterline(waveFolder,frame_seq,head_pos);
DrawCenterline(waveFolder,frame_seq);

% Check if head-tail directions of centerlines are right
prompt = 'Reverse centelines? 0 = skip, other values = yes: ';
reverse_flag = input(prompt);
if reverse_flag ~= 0
    prompt = 'Reverse seq = ?';
    reverse_seq = input(prompt);
    reverse_centerline([waveFolder 'centerline\'],frame_seq(reverse_seq));
end

L = CenterlineLength(Folder,wave_index);
L24Hz = spline((1:length(AVLr_cali))/frame_rateR,L*pixel_size,(1:length(AVLg_cali))/frame_rateG);
LocalizeCenterline;  % Divide the centerline into anterior and posterior halves using a motorneuron as landmark
save([waveFolder 'wormLength.mat'],'L');

%% AVL & DVB tracking

% Pre-segmantation to detect possible neurons
ExtractNeurons(Folder,'green','all')    % set segmentation threshold in 'neuron_seg\NeuronSegConfig.m'

% Tracking performed in GCaMP images
MapCenterlineToOriginSize(Folder,wave_index,rfp_frame_seq);
MapRFPCenterlineToGCaMP(Folder,wave_index);
TrackDVBByCenterline(Folder, gcamp_frame_seq,wave_index,'g');
TrackAVLByLocalSeg(Folder,gcamp_frame_seq,gcamp_frame_seq(1),wave_index,'g');

% Map to the closest wCherry image for neuron positions in wCherry images
MapGCaMP2RFP_single(Folder, 'DVB', wave_index, rfp_frame_seq,gcamp_frame_seq);
MapGCaMP2RFP_single(Folder, 'AVL', wave_index, rfp_frame_seq,gcamp_frame_seq);
 

%% AVL & DVB fluorescence intensities. Use cali if laser field needs to be calibrated
DVBr = SingleNeuronIntensity_RG('DVB','r',[],Folder, wave_index);      % wCherry
DVBr_cali = SingleNeuronIntensity_RG('DVB','r','cali',Folder, wave_index);
DVBg = SingleNeuronIntensity_RG('DVB','g',[],Folder, wave_index);    % GCaMP
DVBg_cali = SingleNeuronIntensity_RG('DVB','g','cali',Folder, wave_index);

AVLr = SingleNeuronIntensity_RG('AVL','r',[],Folder, wave_index);    %wCherry
AVLr_cali = SingleNeuronIntensity_RG('AVL','r','cali',Folder, wave_index);
AVLg = SingleNeuronIntensity_RG('AVL','g',[],Folder, wave_index);    %GCaMP
AVLg_cali = SingleNeuronIntensity_RG('AVL','g','cali',Folder, wave_index);

% normalize using fluorescence intensities in wCherry images
ratioAVL = AVLg_cali./spline((1:length(AVLr_cali))/frame_rateR,AVLr_cali,(1:length(AVLg_cali))/frame_rateG)';
ratioDVB = DVBg_cali./spline((1:length(DVBr_cali))/frame_rateR,DVBr_cali,(1:length(DVBg_cali))/frame_rateG)';

%% Intestinal fluorescence intensities
GreenCenterlineCali(Folder, wave_index,gcamp_frame_seq);
I = WaveIntensity(Folder,wave_index,gcamp_frame_seq,'g');
Ir = WaveIntensity(Folder,wave_index,rfp_frame_seq,'r');
save([waveFolder 'waveIntensity.mat'],'I','Ir');

% smooth and normalize with wCherry
xvalR = (1:length(Ir))/frame_rateR;
xval = (1:length(I))/frame_rateG;

for i = 1:49                        
Ir_24(:,i) = spline(xvalR,Ir(:,i),xval);   
end
sigma = 1;
gausFilter = fspecial('gaussian', [3,3], sigma);
I_ratio = imfilter(I, gausFilter, 'replicate')./imfilter(Ir_24, gausFilter, 'replicate'); % normalize with wCherry
for i = 1:49
I_ratioS(:,i) = smooth(I_ratio(:,i),0.02,'loess');   %smooth
end
for i = 1:49
I_ratioN(:,i) = (I_ratioS(:,i)-min(I_ratioS(:,i)))/range(I_ratioS(:,i)); %normalize to [0,1]
end
save([waveFolder 'waveIntensity.mat'],'I','Ir','Ir_24','I_ratio','I_ratioS','I_ratioN')


