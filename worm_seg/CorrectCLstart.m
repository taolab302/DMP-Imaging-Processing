function CorrectCLstart(OutputFolder)
    load([OutputFolder 'WormRegionPos.mat']);
    centerline_start_xy = load([OutputFolder 'centerline_start_xy.txt']);
    centerline_start = centerline_start_xy(:,[2,1]);
    save([OutputFolder 'WormRegionPos.mat'],'worm_pos','worm_regions','raw_worm_pos','Skip_List','centerline_start');
end