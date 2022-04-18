function anchor = ComputeAVLAnchor(Folder, head_anchor_index, theta_thres, d_thres1, d_thres2, lambda) 
% Compute anchor coordinate of AVL by centerline
PosFolder = [Folder 'neuron_pos\green\'];
NeuronFolder = [Folder 'GCaMP_Neuron\'];
NeuronFiles = dir([NeuronFolder '*.mat']);
FrameSeq = 1:length(NeuronFiles);
anchor = zeros(length(FrameSeq),2);
theta_thres = pi/180 * theta_thres;
head = GetHeadPoints(Folder);
head_anchor = GetQuantilePoints(Folder,head_anchor_index);
HiH1 = head - head_anchor;
for i = 1:length(FrameSeq)
    disp(['Processing ' num2str(i) '/' num2str(length(FrameSeq))]);
    NeuronFile = load([NeuronFolder NeuronFiles(i).name]);
    neurons = NeuronFile.neurons;
    dis = 0; best_index = 0;
    for j = 1:size(neurons,1)
        H1Nx = neurons(j,:) - head(i,:);
        angle = acos(dot(HiH1,H1Nx)/(norm(HiH1)*norm(H1Nx)));
        if angle < theta_thres && norm(H1Nx) > d_thres1 && norm(H1Nx) < d_thres2 && norm(H1Nx) > dis
            dis = norm(H1Nx);
            best_index = j;
        end
    end
    if best_index == 0
        anchor(i,:) = head(i,:);
    else
        anchor(i,:) = neurons(best_index,:);
    end
    anchor(i,:) = (anchor(i,:) + lambda*head(i,:)) / (1 + lambda);
end
anchor(:,[1 2]) = anchor(:,[2 1]);
output_name = [PosFolder 'AVL_anchor.txt'];
fid = fopen(output_name,'w');  
for i = 1:length(FrameSeq)
    fprintf(fid,'%d    %d\n', anchor(i,1), anchor(i,2));
end
fclose(fid);