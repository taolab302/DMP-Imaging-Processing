function MapRFPCenterlineToGCaMP(Folder,wave_index)
waveFolder = [Folder 'Wave\wave-' num2str(wave_index) '\'];
OutFolder = [waveFolder 'centerline_GCaMP\'];
CenterlineFolder = [waveFolder 'centerline_origin\'];
SyncFile = load([Folder 'sync_struc.mat']);
MatchIndex = SyncFile.sync_struc.match_index;
if ~exist(OutFolder, 'dir')
    mkdir(OutFolder);
end
for i = 1:length(MatchIndex)
    disp(['Processing ' num2str(i) '/' num2str(length(MatchIndex))]);
    if exist([CenterlineFolder num2str(MatchIndex(i)) '.mat'],'file')
        load([CenterlineFolder num2str(MatchIndex(i)) '.mat']);
        OutName = [OutFolder num2str(i) '.mat'];
        save(OutName, 'centerline');
    end
end