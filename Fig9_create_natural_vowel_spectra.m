% function Fig9_create_natural_vowel_spectra(saveFig, LatexDir)
function Fig9_create_natural_vowel_spectra(saveFig, LatexDir)

if ~exist('saveFig', 'var')
    saveFig= 0;
end
if ~exist('LatexDir', 'var')
    LatexDir= ['figures' filesep];
end

dirStruct.loading_dir= ['data' filesep 'DanishData' filesep];
dirStruct.Root_savingDir='/media/parida/DATAPART1/Matlab/DropboxOutput/LF_speech_analysis/Output/Danish/';
dirStruct.latexDir= LatexDir;

% chinID=361;

if ~exist('chinID', 'var')
    allfiles=dir([dirStruct.loading_dir '*.mat']);
else
    allfiles=dir([dirStruct.loading_dir '*' num2str(chinID) '*']);
end

allChinSpikeData = [];
for chinVar=1:length(allfiles)
    temp = load([dirStruct.loading_dir allfiles(chinVar).name]);
    allChinSpikeData = [allChinSpikeData; temp.spike_data']; %#ok<AGROW>
end

doPlot= 1;

chin_track_unit_spl= unique([ [allChinSpikeData.chinID]', [allChinSpikeData.track]', [allChinSpikeData.unit]', [allChinSpikeData.SPL]'], 'rows');
tStart = .48; tEnd = .58;


timeSubDir=sprintf('%stStart%.0fms_tEnd%.0fms%s', dirStruct.Root_savingDir, tStart*1e3, tEnd*1e3, filesep);

for unitVar= 144 % Example = 144th unit [Q355, track#5, unit#8]
    curChinID= chin_track_unit_spl(unitVar, 1);
    curTrack= chin_track_unit_spl(unitVar, 2);
    curUnit= chin_track_unit_spl(unitVar, 3);
    dB_SPL= chin_track_unit_spl(unitVar, 4); 
    
    fName=sprintf('%sLF_clean_speech_Q%d_t%d_u%d_tStart%.0fms_tEnd%.0fms_%.0fdBSPL', timeSubDir, curChinID, curTrack, curUnit, tStart*1e3, tEnd*1e3, dB_SPL);
    if ~exist(fName, 'file')
        helper.LF_speech_analysis_danish(curChinID, curTrack, curUnit, dB_SPL, tStart, tEnd, fName, saveFig, dirStruct, doPlot);
    end
end