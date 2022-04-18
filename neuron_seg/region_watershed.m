function L = region_watershed(img)

h = fspecial('sobel');
fd = double(img);

g = sqrt(imfilter(fd, h, 'replicate').^2 + imfilter(fd, h', 'replicate').^2);
figure;imagesc(g);axis image;
L = watershed(g);

end