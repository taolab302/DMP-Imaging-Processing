function ExtractNeurons(Folder,channel,frame_seq)
% Extract neurons from image in Folder

channel = lower(channel);
if strcmp(channel,'g')==1 || strcmp(channel,'green')==1
    ImgFolder = [Folder 'GCaMP\'];
    OutFolder = [Folder 'GCaMP_Neuron\'];
elseif strcmp(channel,'r')==1 || strcmp(channel,'red')==1
    ImgFolder = [Folder 'RFP\'];
    OutFolder = [Folder 'RFP_Neuron\'];
end
if ~exist(OutFolder,'dir')
    mkdir(OutFolder);
end

image_format = '.tiff';
img_seq = GetImageSeq(ImgFolder,image_format);
image_time = img_seq.image_time;
image_prefix = img_seq.image_name_prefix;

if strcmp(frame_seq, 'all') == 1
    frame_seq = 1:length(image_time);
end

tic
for i=1:length(frame_seq)
	image_name = [image_prefix num2str(image_time(frame_seq(i))) image_format];
    disp([num2str(frame_seq(i)) ': ' image_name]);
    
	img = imread([ImgFolder image_name]);

	% Extract neurons in image
	[worm_region,neurons,intensities] = ExtractNeuronInImage(img);

	% Save the neurons data
	output_name = [OutFolder image_name '.mat'];
	save(output_name, 'worm_region', 'neurons','intensities');

% 	% Update neuron intensity threshold
% 	if i~=1
% 		Neuron_Intensity = min(intensities)*0.95;
% 	end
end
time = toc;
disp(['Total time cost: ' num2str(time) ' s']);
end