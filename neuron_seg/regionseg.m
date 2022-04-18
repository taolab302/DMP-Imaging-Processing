function g = regionseg(img,seed_points,dist_thres)
% 利用种子点做区域分割, seed_points为Nx2的数组，表示种子点坐标
% 注：在计算前，需要将种子点合并，即在局部区域中只有一个种子点

seed_num = length(seed_points(:,1));
g = zeros(size(img));
tmp = zeros(size(img));
[height,width] = size(img);
label = 1;
a = 1.5;
br = dist_thres;
region_ratio = 0.2;

for i=1:seed_num
    sr = seed_points(i,1);
    sc = seed_points(i,2);
    min_r = max(1, sr-dist_thres);
    min_c = max(1, sc-dist_thres);
    max_r = min(height, sr+dist_thres);
    max_c = min(width, sc+dist_thres);
    region_img = img(min_r:max_r,min_c:max_c);
    
%     mask = DiskFilter(2*dist_thres+1);
    mask = (fspecial('disk',dist_thres) > 0);
    region_img = region_img.*mask;
    
%     figure;imagesc(region_img);axis image;
    
%     h = fspecial('sobel');
%     fd = double(region);
%     fg = sqrt(imfilter(fd, h, 'replicate').^2 + imfilter(fd, h', 'replicate').^2);
%     figure;imagesc(fg);axis image;
    sr = sr-min_r+1;
    sc = sc-min_c+1;
    seedvalue = region_img(sr,sc);
    small_region = region_img(sr-br:sr+br,sc-br:sc+br);
    nonzero_index = small_region>0;
%     T = max(small_region(nonzero_index)) - ...
%         min(small_region(nonzero_index));
    I_Dif = small_region(nonzero_index) - seedvalue;
    T = mean(abs(I_Dif));
    
    % 如果small region的选择已经覆盖了整个区域
    if T <= region_ratio * seedvalue
        TI = zeros(size(region_img));
        TI(sr-br:sr+br,sc-br:sc+br) = small_region;
        TI( TI > 0) = 1; % non-zero index
    else
        TI = (region_img>0) & (abs(region_img - seedvalue)<=a*T);
    end
     % 保证TI的连通性
    TI = extract_maxconncomp(TI);
    
    tmp(:,:) = 0;
    tmp(min_r:max_r,min_c:max_c) = TI*label;
    g = g+tmp;%将像素叠加到label图像中，而不是block赋值，这样会将原先的值覆盖
    g(g>label) = label;%将重合的部分溶解，认为是最新添加的内容
    label = label+1;
end
end