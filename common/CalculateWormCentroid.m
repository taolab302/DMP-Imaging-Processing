function center = CalculateWormCentroid(img, crop_region)
% crop_region: [min_u, max_y, min_x, max_x]

crop_img = double(img(crop_region(1):crop_region(2),crop_region(3):crop_region(4)));
crop_img = medfilt2(crop_img,[5,5]);
binary_img = crop_img > 150;
binary_img = bwareaopen(binary_img,5,8);
[y,x] = find(binary_img>0);
center = [mean(y) mean(x)] + [crop_region(1) crop_region(3)];

end