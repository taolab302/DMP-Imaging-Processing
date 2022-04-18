function SplitGCaMP_RFP( in_folder, out_folder, mode )
% Separate the images in folder to GCaMP and RFP images
% mode specify the sequence of GCaMP and RFP images

% Make directory for saving GCaMP and RFP images
GCaMP_Folder = [out_folder 'GCaMP'];
RFP_Folder = [out_folder 'RFP'];
if ~exist(GCaMP_Folder,'dir')
    mkdir(GCaMP_Folder);
end
if ~exist(RFP_Folder,'dir')
    mkdir(RFP_Folder);
end

GCAMP_LABEL = 1;
RFP_LABEL = 2;
cycle = length(mode);
image_names = dir([in_folder '*.tiff']);
image_num = length(image_names);
group_num = floor(image_num/cycle);
for i=1:group_num
    % Svaing tne images per group
    disp(['Processing ' num2str(i) '-th group/' num2str(group_num)]);
    
    for j=1:cycle
        index = (i-1)*cycle+j;
        src_imagename = [in_folder char(image_names(index).name)];
        if mode(j) == GCAMP_LABEL
            out_imagename = [GCaMP_Folder '\' char(image_names(index).name)];
            movefile(src_imagename,out_imagename);
        elseif mode(j) == RFP_LABEL
            out_imagename = [RFP_Folder '\' char(image_names(index).name)];
            movefile(src_imagename,out_imagename);
        end
    end
end

% Processing remain images
if length(image_names) > group_num*cycle
    for j=1:(length(image_names)-group_num*cycle)
        index = group_num*cycle + j;
        src_imagename = [in_folder char(image_names(index).name)];
        if mode(j) == GCAMP_LABEL
            out_imagename = [GCaMP_Folder '\' char(image_names(index).name)];
            movefile(src_imagename,out_imagename);
        elseif mode(j) == RFP_LABEL
            out_imagename = [RFP_Folder '\' char(image_names(index).name)];
            movefile(src_imagename,out_imagename);
        end
    end
end
end

