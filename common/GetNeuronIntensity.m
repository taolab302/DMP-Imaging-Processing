function [intensity,ratio] = GetNeuronIntensity(img,neuron_pos,neuron_radius,intensity_ratio)
% Get neuron mean intensity
% Neuorn Pos: [x,y]

neuron_pos = int32(neuron_pos);
[height,width] = size(img);

I = zeros(1,ceil(neuron_radius^2));
count = 0;
nx = 0; ny = 0;
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
intensity = mean(sort_I(1:p));
background = mean(sort_I(end-p+1:end)); % p pixles in the bottom of list
ratio = intensity/background;

end