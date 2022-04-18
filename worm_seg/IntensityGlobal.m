dateFolder = 'F:\Fluorescence\20210714\';
% worm_index = 17;
% mapFolder = [dateFolder 'Intestine-' num2str(worm_index) '-Map\'];
% save_name = [dateFolder 'Intestine-' num2str(worm_index) ' Global_I'];
mapFolder = [dateFolder 'F' num2str(worm_index) '\GCaMP_Map\'];
save_name = [dateFolder 'F' num2str(worm_index) '\Global_I'];

ratio = 1/20;

img_size = 2048;

map_imgs = dir([mapFolder '*.tiff']);
img_num = length(map_imgs);
global_I = zeros(img_num,1);
tic
for i = 1:img_num
   map_img = imread([mapFolder, map_imgs(i).name]);
   map_img = sort(reshape(map_img,img_size^2,1),'descend');
   map_img = map_img(1:floor(ratio*img_size^2));
   global_I(i) = mean(map_img);
   disp([dateFolder(11:18) '-F' num2str(worm_index) ': ' num2str(i) '/' num2str(img_num)]);
%    if mod(i,100)==0
%        disp([dateFolder(11:18) '-F' num2str(worm_index) ': ' num2str(i) '/' num2str(img_num)]);
%    end
end
toc
save(save_name,'global_I');
disp([dateFolder(11:18) '-F' num2str(worm_index) '  Finished!'])