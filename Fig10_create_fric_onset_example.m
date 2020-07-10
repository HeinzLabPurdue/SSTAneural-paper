% function Fig10_create_fric_onset_example(saveFig, LatexDir)
function Fig10_create_fric_onset_example(saveFig, LatexDir)

if ~exist('saveFig', 'var')
    saveFig= 0;
end
if ~exist('LatexDir', 'var')
    LatexDir= ['figures' filesep];
end
dirStruct.latexDir= LatexDir;
dirStruct.loading_Dir= ['data' filesep 'DanishData' filesep];

chinID=355;

if ~exist('chinID', 'var')
    allfiles=dir([dirStruct.loading_Dir '*.mat']);
else
    allfiles=dir([dirStruct.loading_Dir '*' num2str(chinID) '*']);
end

allChinSpikeData = [];
for chinVar=1:length(allfiles)
    temp = load([dirStruct.loading_Dir allfiles(chinVar).name]);
    allChinSpikeData = [allChinSpikeData; temp.spike_data']; %#ok<AGROW>
end

chin_track_unit_spl= unique([ [allChinSpikeData.chinID]', [allChinSpikeData.track]', [allChinSpikeData.unit]', [allChinSpikeData.SPL]'], 'rows');
tStart = .74; tEnd = .85; % Burst window (/s/)

for unitVar=6
    curChinID= chin_track_unit_spl(unitVar, 1);
    curTrack= chin_track_unit_spl(unitVar, 2);
    curUnit= chin_track_unit_spl(unitVar, 3);
    dB_SPL= chin_track_unit_spl(unitVar, 4);
    
    helper.compare_envs_danish(curChinID, curTrack, curUnit, dB_SPL, tStart, tEnd, saveFig, dirStruct);
end
