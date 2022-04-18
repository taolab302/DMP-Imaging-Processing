function NeuronMatch(Folder,frame_seq)

Neuron_Folder = [Folder 'Neuron'];
if ~exist(Neuron_Folder,'dir')
    disp('No Neuron Folder');
    return;
end

Image_Folder = [Folder 'RFP\'];
image_names = dir([Image_Folder '*.tiff']);
image_num = length(frame_seq);

tic
for t=1:image_num
    image_index = frame_seq(t);
    worm_image = imread([Image_Folder image_names(image_index).name]);
    
    neuron_dataname = [Neuron_Folder '\' image_names(image_index).name '.mat'];
    neuron_data = load(neuron_dataname);
    crop_region = neuron_data.crop_region;
    current_worm_center = CalculateWormCentroid(worm_image,crop_region);
    neuron_pos = neuron_data.neuron_pos;
    neuron_num = length(neuron_pos);
    neuron_pos = neuron_pos + repmat([crop_region(1) crop_region(3)],neuron_num,1);
    
    if t==1
        % 计算线虫中心，确定神经元数目
        initail_neuron_num = neuron_num;
        last_worm_center = current_worm_center;
        neurons_tracking_pos = zeros(initail_neuron_num,2,image_num);
        neurons_tracking_pos(:,:,1) = neuron_pos(:,:);
        continue;
    end
    
    worm_offset = current_worm_center - last_worm_center;
    disp([num2str(t) '-th worm offset: ' num2str(worm_offset)]);
    
    for n=1:initail_neuron_num
        last_neuron_pos = neurons_tracking_pos(n,:,t-1);
        last_neuron_pos = last_neuron_pos + worm_offset;
        
        % 选择最近的神经元
        neuron_dist = sum((repmat(last_neuron_pos,neuron_num,1) - neuron_pos).^2,2);
        min_index = find(neuron_dist == min(neuron_dist));
        neurons_tracking_pos(n,:,t) = neuron_pos(min_index(1),:);
    end
    last_worm_center = current_worm_center;
end
time = toc;
disp(['Total time cost: ' num2str(time) ' s']);

% 保存神经元号追踪状态
for n=1:initail_neuron_num
    fid = fopen(['Neuron_' num2str(n)],'wt');  
    for t = 1:image_num
        fprintf(fid,'%d    %d\n',neurons_tracking_pos(n,2,t),neurons_tracking_pos(n,1,t));
    end
    fclose(fid);
end
end