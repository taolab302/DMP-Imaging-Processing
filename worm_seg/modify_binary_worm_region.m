function  b_worm_region = modify_binary_worm_region(binary_worm_region,worm_area)
    worm_region_dist = bwdist(binary_worm_region);
    b_worm_region = worm_region_dist>=4;

    [b_worm_region,~] = Denoise_And_Worm_Locate(b_worm_region, worm_area);  
     b_worm_region = imdilate(b_worm_region,strel('disk',2));
end