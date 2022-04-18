function SyncGCaMP_RFP(img_folder,output_name)
% Synchronize the GCaMP and RFP images

gcamp_seq = GetImageSeq([img_folder 'GCaMP\'], '.tiff');
rfp_seq = GetImageSeq([img_folder 'RFP\'], '.tiff');
sync_struc = SyncImageGroups(gcamp_seq,rfp_seq);
save([img_folder output_name],'sync_struc','gcamp_seq','rfp_seq');

end