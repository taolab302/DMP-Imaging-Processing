function box = boundingbox(bwregion)
% calculate bounding box of binary region

global visited_points;
global point_list;
global point_index;

% Find candidate corner points
% points = detectMinEigenFeatures(bwregion);
% point_locations = points.Location;
point_locations = regionprops(bwregion,'Extrema');
point_locations = point_locations.Extrema;

% Show corner points in image
% imshow(bwregion);hold on;
% plot(point_locations(:,1),point_locations(:,2),'r*');

% Connect points
N = size(point_locations,1);
visited_points = zeros(1,N);
point_list = zeros(1,N);
point_index = 0;

AddPoint(1);
current_point = point_locations(GetCurrentPoint(),:);
candidate_index = find(visited_points == 0);
dist = sum((point_locations(candidate_index,:) - repmat(current_point, length(candidate_index),1)).^2,2);
index = candidate_index(find(dist == min(dist), 1, 'first'));
AddPoint( index );
last_point = current_point;
current_point = point_locations(GetCurrentPoint(),:);
current_dir = current_point - last_point;

alpha = 2;
while (sum(visited_points == 0) > 0)
    candidate_index = find(visited_points == 0);
    point_dirs = point_locations(candidate_index,:) - repmat(current_point, length(candidate_index), 1);
    angles = calculate_angle(point_dirs, current_dir);
    angle_index = angles < pi/2;
    if length(angle_index) == 1
        index = candidate_index(angle_index);
    else
        candidate_index = candidate_index(angle_index);
        dist = sum((point_locations(candidate_index,:) - repmat(current_point, length(candidate_index),1)).^2,2);
        dist = dist + alpha*tan(angles(angle_index));
        index = candidate_index(find(dist == min(dist), 1, 'first'));
    end
    AddPoint( index );
    last_point = current_point;
    current_point = point_locations(GetCurrentPoint(),:);
    current_dir = current_point - last_point;
end
point_list = point_locations(point_list,:);

% Show point list in image
% imshow(bwregion);hold on;
% plot(point_list(:,1),point_list(:,2),'r*-');

% Calculate bounding box
corner_points = FindCornerPoints(point_list, 6);

% Show corner points
box = [corner_points; corner_points(1,:)];
imshow(bwregion);hold on;
plot(box(:,1),box(:,2),'r*-');

end

function AddPoint(index)
    global visited_points;
    global point_list;
    global point_index;

    point_index = point_index + 1;
    visited_points(index) = 1;
    point_list(point_index) = index;
end

function point = GetCurrentPoint()
    global point_list;
    global point_index;
    
    if point_index > 0
        point = point_list(point_index);
    else
        point = [nan,nan];
    end
end

function angles = calculate_angle(point_dirs, base_dir)
    angles = atan2(point_dirs(:,2), point_dirs(:,1)) - atan2(base_dir(2), base_dir(1));
    angles = abs(angles);
    index = angles > pi;
    angles(index) = 2*pi - angles(index);
end

function corner_points = FindCornerPoints(points, n)
    edge1 = [points(end,:) - points(1,:); points(1:end-1,:) - points(2:end,:)];
    edge2 = [points(2:end,:) - points(1:end-1,:); points(1,:) - points(end,:)];
    angles = atan2(edge2(:,2), edge2(:,1)) - atan2(edge1(:,2), edge1(:,1));
    angles = abs(angles);
    index = angles > pi;
    angles(index) = 2*pi - angles(index);
    [~,idx] = sort(angles,'ascend');
    corner_points = points(sort(idx(1:n)),:);
end