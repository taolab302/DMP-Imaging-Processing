function [cc,cc_info] = ExtractAnteriorFrag(binary_img,worm_img, I_thres,frag_center_prev)
%     local_front_low = 450;
%     local_front_high = 10000;
    local_area_low = 300;
    local_area_high = 600;
    
    [height,width] = size(binary_img);
    
    masked_img = binary_img.*worm_img;
    masked_img(~(masked_img>I_thres)) = 0;
    masked_img(masked_img>0) = 1;
%     masked_img = bwmorph(masked_img,'majority');
%     masked_img = imfill(masked_img,'holes');
    cc = bwconncomp(masked_img,4);
    cc_info = zeros(length(cc.PixelIdxList),4);  %col-1:area, col-2:dist to frag_center_prev , col-3&4: frag_center
    for k = 1:length(cc.PixelIdxList)
        fragment = cc.PixelIdxList{k};
        cc_info(k,1) = length(fragment);

        fragment_row = mod(fragment, height);
        fragment_column = ceil(fragment / height);
        pt = [fragment_row,fragment_column];
        if length(fragment)>50                   
            corner_index = convhull(pt);
            corner = pt(corner_index,:);
            frag_center = [mean(corner(1:end-1,1)),mean(corner(1:end-1,2))];
        else
            corner = pt;
            frag_center = [mean(corner(:,1)),mean(corner(:,2))];
        end
        
        frag_center(1) = min(max(frag_center(1),1),height);
        frag_center(2) = min(max(frag_center(2),1),width);
%         frag_center = [mean(corner(1:end-1,1)),mean(corner(1:end-1,2))];
        if ~isempty(frag_center_prev)
            cc_info(k,2) = norm(frag_center - frag_center_prev);
        end
        cc_info(k,3:4) = frag_center;
    end
end