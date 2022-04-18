% Configuration: Centerline
Partition_Num = 49; % Centerline equally divided into 49 segments by default

% Configuration: Binarization
Low_Binary_Thres = 2; % Low threshold for worm region
Hole_Ratio = 0.03;     % Fill the holes within the worm 
BoundaryWidth = 3;    % Width of image edges to be cut off
Frame_Skip_Thres = 0.35; % A frame will be added in Skip list if the area of detected worm is small than 35% of the area in the first frame.

% Configuration: Worm Speed
% Segmentation will be considered wrong if worm displacement between two adjacent frames is larger than MAX_LENGTH.
MAX_LENGTH = 50; % pixels

% Configuration: Omega turn
