function filtered_img = GaussianLowPass(img,sigma)
% Use gaussian low-pass filter inputed image

[height, width] = size(img);
center = ceil(size(img)/2);
[XX, YY] = meshgrid(1:width,1:height);
YY = YY - center(1);
XX = XX - center(2);
low_pass_filter = exp(-(XX.^2 + YY.^2)/(2*sigma^2));
fft_img = fftshift(fft2(img));
filtered_img = abs(ifft2(ifftshift(fft_img .* low_pass_filter)));

end