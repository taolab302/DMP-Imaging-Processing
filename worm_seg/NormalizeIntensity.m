function I_normalize = NormalizeIntensity(I)
    image_num = length(I(:,1));
    seg_num = length(I(1,:));

    sigma = 1;
    gausFilter = fspecial('gaussian', [3,3], sigma);
    I_gauss = imfilter(I, gausFilter, 'replicate');
    I_normalize = zeros(image_num,seg_num);
    for j = 1:seg_num
         % temp = sort(I_gauss(:,i));
         % I_normalize(:,i) = (I_gauss(:,i)-mean(temp(1:25)))/(max(I_gauss(:,i))-mean(temp(1:75)));
         I_gauss(:,j) = smooth(I_gauss(:,j),0.1,'loess');
         I_normalize(:,j) = (I_gauss(:,j)-min(I_gauss(:,j)))/range(I_gauss(:,j));
    end
end