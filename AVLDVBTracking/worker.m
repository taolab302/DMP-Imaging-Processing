Folder = 'H:\Backup\20201210\F3\';
wave_index = 1;
% main_KyIs777;
MapCenterilineToOriginSize(Folder,wave_index,rfp_frame_seq);
MapRFPCenterlineToGCaMP(Folder,wave_index);
TrackDVBByCenterline(Folder, gcamp_frame_seq);
% ComputeAVLAnchor('F:\FluoImages\20201210\F2\','all',3,80,90,210,0.45);

%%