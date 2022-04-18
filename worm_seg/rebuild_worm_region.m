function rebuild_worm = rebuild_worm_region(binary_worm_region,worm_img,slice_width)
    %slice_width = 50;
    step = 15;
    ratio = 0.42;
    [height,width] = size(binary_worm_region);
    slice_numx = floor((width-slice_width)/step)+1;
    slice_numy = floor((height-slice_width)/step)+1;
    rebuild_worm = zeros(height,width);
%     thres = 500;
    area_thres = 100;
%     thres = localmean(worm_img.*~binary_worm_region,1)*2.5;
%     global_thres = median(worm_img.*binary_worm_region,'all');
%     worm_img(worm_img>thres*2) = thres*2;
    thres = localmean(worm_img.*binary_worm_region,0.7);
    % rebuild along x
    for i = 1:slice_numx
        x_start = (i-1)*step+1;
        if i<slice_numx            
            x_end = x_start+slice_width-1;
        else
            x_end = width;
        end
        subimg = worm_img(:,x_start:x_end);
        subworm = binary_worm_region(:,x_start:x_end);
        
        if sum(subworm,'all')>0
            cc = bwconncomp(subworm);
            for k = 1:length(cc.PixelIdxList)
                fragment = cc.PixelIdxList{k};
                fragment_row = mod(fragment, height);
    %             fragment_column = ceil(subworm / height);
                y_start = max(min(fragment_row)-10,1);y_end = min(max(fragment_row)+10,height);
                frag_worm = subworm(y_start:y_end,:);
                local_mean = localmean(frag_worm.*subimg(y_start:y_end,:),0.8);
    %             thres = median(frag_worm.*subimg(y_start:y_end,:),'all');
    %             temp = frag_worm.*subimg(y_start:y_end,:);
    %             local_mean = median(temp(temp>0),'all');
                if local_mean>thres&&length(fragment)>area_thres
                    subrebuild = activecontour(subimg(y_start:y_end,:),frag_worm,20);
                elseif length(fragment)<area_thres
                    subrebuild = zeros(y_end-y_start+1,x_end-x_start+1);
                else
                    subrebuild = frag_worm;
                end
                if sum(subrebuild,'all')>100
                rebuild_worm(y_start:y_end,x_start:x_end) = rebuild_worm(y_start:y_end,x_start:x_end)+subrebuild;
                end
            end     
        end
    end
    
    %rebuild along y
    for i = 1:slice_numy
        y_start = (i-1)*step+1;
        if i<slice_numy            
            y_end = y_start+slice_width-1;
        else
            y_end = height;
        end
        subimg = worm_img(y_start:y_end,:);
        subworm = binary_worm_region(y_start:y_end,:);
        
        if sum(subworm,'all')>0
            cc = bwconncomp(subworm);
            for k = 1:length(cc.PixelIdxList)
                fragment = cc.PixelIdxList{k};
    %             fragment_row = mod(fragment, height);
                fragment_column = ceil(fragment / length(subworm(:,1)));
                x_start = max(min(fragment_column)-10,1);x_end = min(max(fragment_column)+10,width);
                frag_worm = subworm(:,x_start:x_end);
                local_mean = localmean(frag_worm.*subimg(:,x_start:x_end),0.8);
    %             thres = median(frag_worm.*subimg(:,x_start:x_end),'all');
    %             temp = frag_worm.*subimg(:,x_start:x_end);
    %             local_mean = median(temp(temp>0),'all');
                if local_mean>thres||length(fragment)<area_thres
                    subrebuild = activecontour(subimg(:,x_start:x_end),frag_worm,20);
                elseif length(fragment)<area_thres
                    subrebuild = zeros(y_end-y_start+1,x_end-x_start+1);
                else
                    subrebuild = frag_worm;
                end

                if sum(subrebuild,'all')>100
                rebuild_worm(y_start:y_end,x_start:x_end) = rebuild_worm(y_start:y_end,x_start:x_end)+subrebuild;
                end
            end   
        end
    end
       
    rebuild_worm(rebuild_worm>0) = 1;

end

function localmean = localmean(I,ratio)
    temp = sort(I(:));
    temp = temp(temp>0);
    localmean = mean(temp(1:floor(length(temp)*ratio)));
end