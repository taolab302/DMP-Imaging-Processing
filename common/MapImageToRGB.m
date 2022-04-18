function map_image = MapImageToRGB(image,map_range,type)
% Map image into 8-bit
% 
% Input parameters:
% image:     original image
% map_range: data range when making video. If the type of data is UINT8,
%            this value will be [0,255], and if the type of data is UINT16,
%            this value is needed to be set by image, such as [min_v, max_v]

[height, width] = size(image);
min_value = map_range(1);
max_value = map_range(2);
map_image = zeros(height, width, 3);

image = double(image);
image(image>max_value) = max_value;
image(image<min_value) = min_value;
image = uint8((image-min_value)*255/(max_value-min_value));

if strcmp(type,'red') == 1
    map_image(:,:,1) = image;
elseif strcmp(type,'green') == 1
    map_image(:,:,2) = image;
end
map_image = uint8(map_image);

end