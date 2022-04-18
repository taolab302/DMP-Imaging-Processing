function [corrected_img,Background] = CalibrateImage(img,pattern,MaxBackground)
% using image pattern (inverse_GCaMP/inverse_RFP) to modify the original image
% Hint: make sure the object is larger than background!
% 
% For GCaMP image MaxBacground can be 200, for RFP image this value can be 150

BorderLen = 50;
% ratio = 1.1;
ratio = 1;
[height,width] = size(img);
half_height = floor((height+1)/2); half_width = floor((width+1)/2);
img = medfilt2(img,[5,5]);
% figure;imagesc(img);

% search four regions of image to calculate a proper background
backgrounds = zeros(1,4);
corner_img1 = img(1:BorderLen,(half_width-BorderLen+1):half_width);
backgrounds(1) = mean(corner_img1(:));

corner_img2 = img((half_height-BorderLen+1):half_height,(width-BorderLen+1):width);
backgrounds(2) = mean(corner_img2(:));

corner_img3 = img((height-BorderLen+1):height,(half_width-BorderLen+1):half_width);
backgrounds(3) = mean(corner_img3(:));

corner_img4 = img((half_height-BorderLen+1):half_height,1:BorderLen);
backgrounds(4) = mean(corner_img4(:));

for i=1:4
    if backgrounds(i)*ratio > MaxBackground
        backgrounds(i) = nan;
    end
end
backgrounds = backgrounds(~isnan(backgrounds));
if isempty(backgrounds)
    Background = MaxBackground/ratio;
else
    Background = ratio*max(backgrounds);
end

% use pattern to calibrate image
corrected_img = (double(img) - Background) .* (img > Background);
% corrected_img  = corrected_img .* pattern + Background;
corrected_img  = corrected_img .* pattern + Background;

% % sigma = 1;
% % gausFilter = fspecial('gaussian', [3,3], sigma);
% % img_gauss = imfilter(img, gausFilter, 'replicate');
% corrected_img  = double(img).* pattern;

end

