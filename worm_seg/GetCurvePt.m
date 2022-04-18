function curve = GetCurvePt(BW,start_pt)
    pt_num = length(find(BW>0));
    curve = zeros(pt_num,2);
    offsets = [0 -1; -1 -1; -1 0; -1 1; 0 1; 1 1; 1 0; 1 -1]; %clockwise
    curve(1,:) = start_pt;
    direct_index = 4;
%     tail_end = end_pt(end_dist == max(end_dist),:);
    for i = 2:length(curve(:,1))
        if curve(i-1,1)>0
%             if i == 114
%                 pause;
%             end
             for j = direct_index+5:direct_index+12
                 if mod(j,8) == 0
                     dir_num = 8;
                 else
                     dir_num = mod(j,8);
                 end
                 new_pt = curve(i-1,:)+offsets(dir_num,:);
                 if BW(new_pt(1),new_pt(2))>0&&isempty(find(curve(:,1)==new_pt(1)&curve(:,2)==new_pt(2), 1))
                     curve(i,:) = new_pt;
                     direct_index = mod(j,8);
                    break;
                 end   
             end
%              disp(num2str(i))
         end
    end
    curve = curve(curve(:,1)>0,:);
end