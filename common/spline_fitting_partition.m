function [interpolated_line] = spline_fitting_partition(before_interpolate, PARTITION_NUM)
% 分割中心线，获取中心点坐标(interpolated_line)
% 参数说明：before_interpolate 表示获取的边界坐标，PARTITION_NUM 表示分割的数目

POINT_NUM = size(before_interpolate,1);
border_index_m = before_interpolate(:,1);
border_index_n = before_interpolate(:,2);
spline_arc = zeros(POINT_NUM,1);
interpolated_line = zeros(PARTITION_NUM+1,2);

pp_m = csaps(1 : POINT_NUM, border_index_m);                            % doing the cubic smoothing spline for border x index
pp_n = csaps(1 : POINT_NUM, border_index_n);                            % doing the cubic smoothing spline for border y index
curve_m = fnval(pp_m, 1:POINT_NUM);                                     % giving the x location of spline curve
curve_n = fnval(pp_n, 1:POINT_NUM);                                     % giving the y location of spline curve

curve_m_diff = curve_m(2:POINT_NUM) - curve_m(1:POINT_NUM-1);
curve_n_diff = curve_n(2:POINT_NUM) - curve_n(1:POINT_NUM-1);
curve_diff = sqrt(curve_m_diff.^2 + curve_n_diff.^2);
spline_length = sum(curve_diff);

spline_arc(2:POINT_NUM) = cumsum(curve_diff);
spline_pp_m = csaps(spline_arc,curve_m);
spline_pp_n = csaps(spline_arc,curve_n);

interpolated_line(:,1) = fnval(spline_pp_m, 0:spline_length / PARTITION_NUM:spline_length);
interpolated_line(:,2) = fnval(spline_pp_n, 0:spline_length / PARTITION_NUM:spline_length);
end