function filtered_img = GaussianHighPass(img,sigma)
% Use gaussian high-pass filter inputed image

[height, width] = size(img);
center = ceil(size(img)/2);
[XX, YY] = meshgrid(1:width,1:height);
YY = YY - center(1);
XX = XX - center(2);
high_pass_filter = 1 - exp(-(XX.^2 + YY.^2)/(2*sigma^2));

fft_img = fftshift(fft2(img));
filtered_img = abs(ifft2(ifftshift(fft_img .* high_pass_filter)));

end

