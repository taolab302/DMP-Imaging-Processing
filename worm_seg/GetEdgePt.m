function [left,right] = GetEdgePt(edge,dir,centerline_pt)
    %沿法线向两边延伸获得左右轮廓，正方向是右边缘(dorsal)，负方向是左边缘（ventral)
    % dir 中心线方向
    boundary_image = imfill(edge,'holes');
    width = 1000;
    dir = dir/norm(dir);
    dir_p = [-dir(2) dir(1)];
    delta_x = 2*(dir(1)>0)-1;
    delta_y = 1-2*(dir(2)>0);
     %右边缘点
    for j=0:width        
        point = j*dir_p + centerline_pt;
        p = int32(point);
        if boundary_image(p(1),p(2))==1 || boundary_image(p(1)+delta_x,p(2))==1 ||...
           boundary_image(p(1),p(2)+delta_y)==1 || boundary_image(p(1)+delta_x,p(2)+delta_y)==1 ||...
           boundary_image(p(1)-delta_x,p(2))==1 || boundary_image(p(1)-delta_x,p(2)+delta_y)==1
            right = int32(point + [0.5,0.5].*[delta_x,delta_y]);
            break;
        end       
    end
    %左边缘点
    
    delta_x = -delta_y;
    delta_y = delta_x;
    for j=0:width
        point = j*(-dir_p) + centerline_pt;
        p = int32(point);
        if boundary_image(p(1),p(2))==1 || boundary_image(p(1)+delta_x,p(2))==1 ||...
           boundary_image(p(1),p(2)+delta_y)==1 || boundary_image(p(1)+delta_x,p(2)+delta_y)==1 ||...
           boundary_image(p(1)-delta_x,p(2))==1 || boundary_image(p(1)-delta_x,p(2)+delta_y)==1
            left = int32(point + [0.5,0.5].*[delta_x,delta_y]);
            break;
        end
    end
end


end