function pos = TrackAVLByAnchor(Folder, FrameSeq) 
% Track AVL by Anchor
PosFolder = [Folder 'neuron_pos\green\'];
NeuronFolder = [Folder 'GCaMP_Neuron\'];
NeuronFiles = dir([NeuronFolder '*.mat']);
anchor_pos = load([PosFolder 'AVL_anchor.txt']);
if strcmp(FrameSeq, 'all') == 1
    FrameSeq = 1:length(NeuronFiles);
end
pos = zeros(length(FrameSeq),2);
for ii = 1:length(FrameSeq)
    i = FrameSeq(ii);
    if length(FrameSeq) == 1
        disp('Replace neuron position by extraction.');
    else
        disp(['Processing ' num2str(ii) '/' num2str(length(FrameSeq))]);
    end
    NeuronFile = load([NeuronFolder NeuronFiles(i).name]);
    neurons = NeuronFile.neurons;
    dis = inf; best_index = 0;
    for j = 1:size(neurons,1)
        neuron = [neurons(j,2) neurons(j,1)];
        nowdis = dist(anchor_pos(i,:),neuron');
        if nowdis < dis
            dis = nowdis;
            best_index = j;
        end
    end
    if best_index == 0
        pos(ii,1) = anchor_pos(i,2);
        pos(ii,2) = anchor_pos(i,1);
    else
        pos(ii,1) = neurons(best_index,2);
        pos(ii,2) = neurons(best_index,1);
    end
end

if length(FrameSeq) > 1
    output_name = [PosFolder 'AVL.txt'];
    fid = fopen(output_name,'w');  
    for i = 1:length(FrameSeq)
        fprintf(fid,'%d    %d\n', pos(i,1), pos(i,2));
    end
    fclose(fid);
end

    



        