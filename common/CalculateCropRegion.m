function region = CalculateCropRegion(img,p)
% Crop image that only contains worm
% Input
% p: parameters in croping
% region: [y_min, y_max, x_min, x_max]

margin = 20;

% Binarize image
binary_worm = img>p;

% Obtain worm region
[height,width] = size(img);
[worm_row, worm_column] = find(binary_worm>0);
y_min = max(min(worm_row) - margin, 1);
y_max = min(max(worm_row) + margin, height);
x_min = max(min(worm_column) - margin, 1);
x_max = min(max(worm_column) + margin, width);
region = [y_min, y_max, x_min, x_max];

end