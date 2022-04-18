function image_seq = GetImageSeq(Image_Dir,image_format)
% Read image filenames

image_filenames = dir([Image_Dir '*' image_format]);
image_num = length(image_filenames);

image_name_prefix = '';
time = zeros(image_num,1);
for i=1:image_num
    image_name = image_filenames(i).name;
    
    % processing char array of image name
    split_pattern = strfind(image_name,'_');
    if isempty(split_pattern)
       split_pattern = 0;
    else
        split_pattern = split_pattern(end);
    end
    str_time = image_name(split_pattern+1:(length(image_name)-length(image_format)));
    time(i) = str2double(str_time);
        
    if i==1
        image_name_prefix = image_name(1:split_pattern);
    end
end
% sort time
image_time = sort(time);

image_seq.image_time = image_time;
image_seq.image_name_prefix = image_name_prefix;
end

