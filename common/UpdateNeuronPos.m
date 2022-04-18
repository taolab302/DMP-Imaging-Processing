function [xnext,ynext] = UpdateNeuronPos(xprev,yprev,search_interval,Intensity_Ratio,Wimage)
% Update the neuron position by its centroid

Cx = zeros(1,ceil(search_interval*search_interval));
Cy = zeros(1,ceil(search_interval*search_interval));
Ci = zeros(1,ceil(search_interval*search_interval));

k = 0;
for x=(xprev-search_interval):(xprev+search_interval)
    for y=(yprev-search_interval):(yprev+search_interval)
         if (x-xprev)^2+(y-yprev)^2 <= search_interval^2
            k = k+1;
            Cx(k) = min(max(1, int32(x)),length(Wimage(1,:)));
            Cy(k) = min(max(1, int32(y)),length(Wimage(:,1)));
            Ci(k) = Wimage(Cy(k),Cx(k));
         end
    end
end
Pixels_Num = k; % number of pixels in search region (circle with radius search interval)
expected_pixel_num = ceil(Pixels_Num*Intensity_Ratio);
[sort_Ci, IDx] = sort(Ci(1:Pixels_Num),'descend');
UpCx = Cx(IDx(1:expected_pixel_num));
UpCy = Cy(IDx(1:expected_pixel_num));
UpCi = sort_Ci(1:expected_pixel_num);

% Update center of mass
UpIc_Energy = sum(UpCi);
xnext = sum(UpCx.*UpCi)/UpIc_Energy;
ynext = sum(UpCy.*UpCi)/UpIc_Energy;
% disp(['Center_x = ',num2str(Center_x),' Center_y = ',num2str(Center_y)]);
end