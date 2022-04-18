function [ca,ba] = ExtractFluoEnergyAndBackground(img,neuron_pos,neuron_radius,intensity_ratio)
% Neuorn Pos: [x,y]

[height,width] = size(img);
I = zeros(1,ceil(neuron_radius^2));

count = 0;
for x = (neuron_pos(1)-neuron_radius):(neuron_pos(1)+neuron_radius)
    for y = (neuron_pos(2)-neuron_radius):(neuron_pos(2)+neuron_radius)
        if ((x-neuron_pos(1))^2+(y-neuron_pos(2))^2) <= neuron_radius^2
            nx = max(1, min(int32(x), width));
            ny = max(1, min(int32(y), height));
            count = count+1;
            I(count) = img(ny,nx);
        end
    end
end
I = I(1:count);

sort_I = sort(I,'descend');
p = int32(count*intensity_ratio);
ca = mean(sort_I(1:p));% intensity threshold for calculating neuron activity
ba = mean(sort_I(p+1:end));
% ca = sum(sum(I));
% ba = 0; 

end