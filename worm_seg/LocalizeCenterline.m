% waveFolder = 'F:\Fluorescence\20201210\F3\Wave\wave-1\';
channel = 'red';
neuron_pos = load([waveFolder 'neuron_pos\' channel '\neuron 03.txt']);  %load position of the landmark
frame_seq = rfp_frame_seq;
% frame_seq = 1617:2001;
% waveFolder = [Folder 'wave\Wave-' num2sr(i) '\'];
L_anterior = zeros(length(frame_seq),1);
L_posterior = zeros(length(frame_seq),1);
for i = frame_seq
%     load([waveFolder 'centerline_GCaMP_cali\' num2str(i) '.mat']);
    load([waveFolder 'centerline_origin\' num2str(i) '.mat']);
    dist_cl = sqrt(sum((neuron_pos(i-frame_seq(1)+1,[2,1])-centerline).^2,2));
    [~,id] = sort(dist_cl);
    if id(1)<id(2)
        c1 = id(1);
        c2 = id(2);        
    else
        c1 = id(2);
        c2 = id(1);
    end
    seg_cl = dot(centerline(c2,:)-centerline(c1,:), neuron_pos(i-frame_seq(1)+1,[2,1])-centerline(c1,:))/norm(centerline(c2,:)-centerline(c1,:));
    L_anterior(i-frame_seq(1)+1) = sum(sqrt(sum((centerline(2:c1,:)-centerline(1:c1-1,:)).^2,2)))+seg_cl;
    L_posterior(i-frame_seq(1)+1) = sum(sqrt(sum((centerline(c1+1:end,:)-centerline(c1:end-1,:)).^2,2)))-seg_cl;
    
end

save([waveFolder 'wormLength.mat'],'L','L24Hz','L_anterior','L_posterior')
figure;plot(L_anterior)
hold on;plot(L_posterior)