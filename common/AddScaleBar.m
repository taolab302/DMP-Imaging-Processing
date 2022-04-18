function img = AddScaleBar(img,scalebar_size)
% Add scale bar in image
% 

height = size(img,1);
width = size(img,2);

if size(img,3) < 3
	new_img = zeros(height, width, 3);
	new_img(:,:,1) = img;
	new_img(:,:,2) = img;
	new_img(:,:,3) = img;
	img = new_img;
end

row_start = floor(0.918*height);
row_end = row_start + 6;
col_start = floor(0.8481*width);
col_end = col_start + scalebar_size-1;

img(row_start:row_end,col_start:col_end,1) = 255;
img(row_start:row_end,col_start:col_end,2) = 255;
img(row_start:row_end,col_start:col_end,3) = 255;

imshow(img);
end