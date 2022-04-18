function ReplaceNeuronPos(Folder,tracking_index, insert_pos, FlouType)
% Split Map text and append into neuron position files

posFolder = [Folder 'neuron_pos\' FlouType '\'];
if strcmp(FlouType,'red')
    neurons_pos = load([posFolder,'RFP_Map.txt']);
elseif strcmp(FlouType,'green')
    neurons_pos = load([posFolder,'GCaMP_Map.txt']);
end
AVLDVBflag = strcmp(tracking_index, 'AVL') || strcmp(tracking_index, 'DVB');
if AVLDVBflag
    neuron_num = 1;
else
    neuron_num = length(tracking_index);
end
Tracking_Length = length(neurons_pos(:,1))/neuron_num;
insert_pos = insert_pos(1);%make sure only insert data in one position

for i = 1:neuron_num
    if AVLDVBflag
        output_name = [posFolder tracking_index '.txt'];
    else
        output_name = [posFolder sprintf('neuron %02d',tracking_index(i)) '.txt'];
    end
    original_neuron_pos = load(output_name);
    updated_neuron_pos = zeros(length(original_neuron_pos), 2);

    % insert neuron positions
    if ~isempty(insert_pos)
        updated_neuron_pos(1:insert_pos-1,:) = original_neuron_pos(1:insert_pos-1,:);
        for k=1:Tracking_Length
            updated_neuron_pos(insert_pos+k-1,:) = neurons_pos((k-1)*neuron_num+i, :);
        end
        % replace neuron postions
        updated_neuron_pos((insert_pos+Tracking_Length):end,:) = original_neuron_pos((insert_pos+Tracking_Length):end,:);
    else
        % append the neuron position
        insert_pos = length(original_neuron_pos) + 1;
        updated_neuron_pos(1:(insert_pos-1),:) = original_neuron_pos(:,:);
        for k=1:Tracking_Length
            updated_neuron_pos(insert_pos+k-1,:) = neurons_pos((k-1)*neuron_num+i, :);
        end
    end

    % write neuron positions into file
    if AVLDVBflag
        write_neuronpos(posFolder,tracking_index,updated_neuron_pos(:,1), updated_neuron_pos(:,2),'w');
        disp(['Succeed: ' tracking_index]);
    else
        write_neuronpos(posFolder,tracking_index(i),updated_neuron_pos(:,1), updated_neuron_pos(:,2),'w');
        disp(['Succeed:  neuron ',num2str(tracking_index(i))]);
    end
end

end