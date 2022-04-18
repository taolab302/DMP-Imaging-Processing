function [img_amend,rest_I] = AmendCenterI(img,height, width,center)
    amend_ratio = 0.3;
    img_size = length(img(:,1));
    radium = ceil(img_size*amend_ratio);
    theta = linspace(0,2*pi,50);
    
    Mask = poly2mask(radium*cos(theta)+center(1),radium*sin(theta)+center(2),height,width);
    Mask_pt = sum(sum(Mask));
    rest_pt = height*width - Mask_pt;
    Masked_img = Mask.*img;
    unMasked_img = (~Mask).*img;
    Amend_I = sort(reshape(Masked_img(Masked_img>0),1,Mask_pt));
%     rest_I = sort(reshape(unMasked_img(unMasked_img>0),1,rest_pt));
    rest_I = sort(reshape(unMasked_img(unMasked_img>0),1,length(unMasked_img(unMasked_img>0))));
    
    Amend_I = median(Amend_I);
    rest_I = median(rest_I);
%     Amend_I = mean(Amend_I(1:floor(Mask_pt*amend_ratio)));
%     rest_I = mean(rest_I(1:floor(rest_pt*amend_ratio)));
    
    delta_I = mean(Amend_I - rest_I);
    img_amend = img - delta_I*Mask.*(img>rest_I);

end