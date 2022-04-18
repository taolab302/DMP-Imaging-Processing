function centroid = CalculateWormCentroid(img)
% Calculate worm centroid

[height, width] = size(img);
[XX,YY] = meshgrid(1:width,1:height);
energy = sum(sum(img));
centroid_x = sum(sum(img.*XX))/energy;
centroid_y = sum(sum(img.*YY))/energy;
centroid = [centroid_y,centroid_x];

end