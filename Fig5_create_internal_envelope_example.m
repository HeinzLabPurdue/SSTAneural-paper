% function Fig5_create_internal_envelope_example(saveFig, LatexDir)
function Fig5_create_internal_envelope_example(saveFig, LatexDir)

if ~exist('saveFig', 'var')
    saveFig= 0;
end
if ~exist('LatexDir', 'var')
    LatexDir= ['figures' filesep];
end
dirStruct.latexDir= LatexDir;
dirStruct.loading_Dir=['data' filesep 'EnglishData' filesep];
dirStruct.OutFig_Dir= '/media/parida/DATAPART1/Matlab/DropboxOutput/hilbert_envelopes/English/';
allfiles=dir([dirStruct.loading_Dir '*.mat']);

allChinSpikeData = [];
for chinVar=1:length(allfiles)
    temp = load([dirStruct.loading_Dir allfiles(chinVar).name]);
    allChinSpikeData = [allChinSpikeData; temp.spike_data']; %#ok<AGROW>
end

all_chin_track_unit_spl= [ [allChinSpikeData.ChinID]', [allChinSpikeData.track]', [allChinSpikeData.unit]', [allChinSpikeData.SPL]'];
chin_track_unit_spl= unique(all_chin_track_unit_spl, 'rows');
tStart3 = 0; tEnd3 = 2.28;


anl.fs= 2e3;
anl.binWidth= 1/anl.fs;
anl.BP.freqs= [2 4 8 16 32 64 128];
anl.BP.order= 2;
anl.stimFs= 20e3;
anl.LP.co_hz= 32;
anl.LP.order= 2;

filt_obj_LP= helper.get_filter_designfilt('lp', anl.LP.co_hz, anl.fs, anl.LP.order);
filt_abj_BP= cell(length(anl.BP.freqs), 1);
for fmVar=1:length(anl.BP.freqs)
    cur_fm= anl.BP.freqs(fmVar);
    filt_abj_BP{fmVar}= helper.get_filter_designfilt('bp', cur_fm*[1/sqrt(2) sqrt(2)], anl.fs, anl.BP.order, 0);
end


for unitVar= 16
    curChinID= chin_track_unit_spl(unitVar, 1);
    curTrack= chin_track_unit_spl(unitVar, 2); %#ok<*PFBNS>
    curUnit= chin_track_unit_spl(unitVar, 3);
    dB_SPL= chin_track_unit_spl(unitVar, 4);
    
    unit_inds= find(ismember(all_chin_track_unit_spl, chin_track_unit_spl(unitVar, :),'rows'));
    if ~isempty(unit_inds)
        curUnitData= allChinSpikeData(unit_inds(1));
        
        curSpikes_pos=cell2mat(curUnitData.SpikeTrains{1,1}');
        curSpikes_neg=cell2mat(curUnitData.SpikeTrains{1,2}');
        
        fName=sprintf('%sHilbEnv_clean_speech_Q%d_t%d_u%d_%.0fdBSPL.png', dirStruct.OutFig_Dir, curChinID, curTrack, curUnit, dB_SPL);
        helper.create_intrinsic_envelopes(curSpikes_pos, curSpikes_neg, anl, filt_obj_LP, filt_abj_BP, fName, dirStruct, saveFig);
    end
end
