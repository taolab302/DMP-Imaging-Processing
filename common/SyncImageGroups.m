function sync_struc = SyncImageGroups(images_seq1,images_seq2)
% Synchronize the GCaMP and RFP images
% 
% Input parameters:
% images_seq1: the sequence of image group 1
% images_seq2: the sequence of image group 2
% Hint: the number of images_seq1 is more than images_seq2
%
% Output parameters:
% sync_struc: synchronous result of image group-1 and group-2

close all;
image_format = '.tiff';
group1_prefix = images_seq1.image_name_prefix;
group1_time = images_seq1.image_time;
group2_prefix = images_seq2.image_name_prefix;
group2_time = images_seq2.image_time;

% Plot interval of two continunous images in each group
% figure;plot(diff(group1_time),'b.');title('Group-1 image time');ylabel('ms');
% figure;plot(diff(group2_time),'b.');title('Group-2 image time');ylabel('ms');

Num = length(group1_time); % #{group1_time} > #{group2_time}
group1_interval = mean(diff(group1_time));
group2_interval = mean(diff(group2_time));
sync_names = cell(Num,2); %The first column is group-1 and the second is group-2
match_index = zeros(Num,1);

Time_Thres = group1_interval/5;
for i=1:Num
    time_1 = group1_time(i);
    sync_names{i,1} = [group1_prefix num2str(time_1) image_format];
    
    time_diff = group2_time - time_1;
    match_flag = abs(time_diff) < (group1_interval + Time_Thres);
    index = find(match_flag == 1,1); 
    if isempty(index)
        sync_names{i,2} = '';
        match_index(i) = nan;
    else
        time_2 = group2_time(index);
        sync_names{i,2} = [group2_prefix num2str(time_2) image_format];
        match_index(i) = index;
    end 
end

% Fill the empty one in sync_names
if isnan(match_index(1))
    disp('The first image cannot be matched, please check!');
    disp('Now we assume the first images in each group are matched.');
    match_index(1) = 1;
    sync_names{1,2} = [group2_prefix num2str(group2_time(1)) image_format];
end

for i=1:Num
    if i==1
        last_match_name = sync_names(1,2);
        last_index = match_index(1);
    elseif ~isempty(char(sync_names(i-1,2)))
        last_match_name = sync_names(i-1,2);
        last_index = match_index(i-1);
    end

    if isempty(char(sync_names(i,2)))
        sync_names(i,2) = last_match_name;
        match_index(i) = last_index;
    end
end

sync_struc.sync_names = sync_names;
sync_struc.match_index = match_index;
sync_struc.interval1 = group1_interval;
sync_struc.interval2 = group2_interval;
end