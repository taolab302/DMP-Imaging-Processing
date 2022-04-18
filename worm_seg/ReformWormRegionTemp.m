function ReformWormRegionTemp(waveFolder)
    orig_format = '.tiff.tif';
    aim_format = '.tiff';
    imgs = dir([waveFolder 'worm_region\*',orig_format]);
    img_num = length(imgs);
    for i = 1:img_num
        orig_name = [waveFolder 'worm_region\' imgs(i).name];
        aim_name = [orig_name(1:(end-length(orig_format))),aim_format];
        copyfile(orig_name, aim_name);
        delete(orig_name);
    end

end