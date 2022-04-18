function h = CreateGaborFilter(w,theta)
% Create Gabor filter
%
% w: wavelength in Gabor filter
% theta: orientation in gabor filter

g = gabor(w,theta);
hsize = floor(1.5*w);
h_full = real(g.SpatialKernel);
cy = (size(h_full,1)+1)/2;
cx = (size(h_full,2)+1)/2;
h = h_full(cy-hsize:cy+hsize,cx-hsize:cx+hsize);
% h = h_full;
end