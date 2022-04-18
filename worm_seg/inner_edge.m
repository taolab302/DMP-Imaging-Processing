function boundary_image = inner_edge(bw_img)

img_dist = bwdist(~bw_img);
boundary_image = (img_dist <= 3 & img_dist > 1.5);

end