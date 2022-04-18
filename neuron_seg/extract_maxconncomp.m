function out = extract_maxconncomp(im)
% ��ȡ�����ͨ��֧

conncomps = bwconncomp(im,8);
pixels_num = zeros(conncomps.NumObjects,1);
out = im;
for i=1:conncomps.NumObjects
    pixels_num(i) = length(conncomps.PixelIdxList{i});
end
max_index = find(pixels_num == max(pixels_num),1);
for i=1:conncomps.NumObjects
    if i~=max_index
        out(conncomps.PixelIdxList{i}) = 0;        
    end
end
% �������ͨ��֧����hole�������
out = imfill(out,'holes');
end
