function [newpos,neuron_intensity,background] = UpdateNeuronData(prev_pos,search_interval,intensity_ratio,Wimage)
% Update the neuron position, instensity and background
% Input parameters:
% prev_pos: [x,y] neuron position
% search_interval: search range
% Intensity_Ratio: the ratio of foreground(neuron) pixels number in the search region
% Wimage: worm image
%
% Output parameters:
% newpos: updated neuron position, [x,y] format
% neuron_intensity: mean neuron intensity
% background: mean intensity of background

Cx = zeros(1,ceil(search_interval*search_interval));
Cy = zeros(1,ceil(search_interval*search_interval));
Ci = zeros(1,ceil(search_interval*search_interval));

k = 0;
xprev = prev_pos(1); yprev = prev_pos(2);
for x=(xprev-search_interval):(xprev+search_interval)
    for y=(yprev-search_interval):(yprev+search_interval)
         if (x-xprev)^2+(y-yprev)^2 <= search_interval^2
            k = k+1;
            Cx(k) = max(1, int32(x));
            Cy(k) = max(1, int32(y));
            Ci(k) = Wimage(Cy(k),Cx(k));
         end
    end
end
Pixels_Num = k; % number of pixels in search region (circle with radius search interval)
expected_pixel_num = ceil(Pixels_Num*intensity_ratio);
[sort_Ci, IDx] = sort(Ci(1:Pixels_Num),'descend');
UpCx = Cx(IDx(1:expected_pixel_num));
UpCy = Cy(IDx(1:expected_pixel_num));
UpCi = sort_Ci(1:expected_pixel_num);

% Update center of mass
newpos = prev_pos;
UpIc_Energy = sum(UpCi);
newpos(1) = sum(UpCx.*UpCi)/UpIc_Energy;
newpos(2) = sum(UpCy.*UpCi)/UpIc_Energy;
neuron_intensity = mean(UpCi);
background = mean(sort_Ci(expected_pixel_num:end));
end