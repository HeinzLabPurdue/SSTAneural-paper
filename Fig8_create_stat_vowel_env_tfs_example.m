% function Fig8_create_stat_vowel_env_tfs_example(saveFig, LatexDir)
function Fig8_create_stat_vowel_env_tfs_example(saveFig, LatexDir)

if ~exist('saveFig', 'var')
    saveFig= 0;
end 
if ~exist('LatexDir', 'var')
    LatexDir= ['figures' filesep];
end
DirStruct.latexDir= LatexDir;
DirStruct.INdata= ['data' filesep];

anl.AcousticDelay= 8e-3; % Approximate stimulus-response delay. 4.59 ms for inv-FIR filtering. ~3 ms for default circuit.
anl.stimDuration= 188e-3;
anl.stimPeriod= 350e-3;
anl.fs= 10e3;
anl.nw= 1.5; % time-bandwidth half product for pmtm
anl.tEdge_hist= (0:1/anl.fs:anl.stimPeriod)+anl.AcousticDelay;
tPlot= (anl.tEdge_hist(1:end-1) + anl.tEdge_hist(2:end))/2-anl.AcousticDelay;
anl.tStart= 38e-3; %38e-3;
anl.tEnd= 188e-3; %anl.stimDuration; % anl.stimDuration || (38e-3+64e-3);

anl.tMask= tPlot>anl.tStart & tPlot<anl.tEnd;
anl.feature_to_use_f0= 'RAW';

% Load saved data

ChinID= 374;
SynCapData= load(sprintf('%sQ%d_allSyncCap.mat', DirStruct.INdata, ChinID));
cur_unit_data=SynCapData.SynCapData(4);


uRatePos= histcounts(cell2mat(cur_unit_data.stat.RAW.pos'), anl.tEdge_hist);
uRateNeg= histcounts(cell2mat(cur_unit_data.stat.RAW.neg'), anl.tEdge_hist);

[sig_raw, fsOrg]= audioread(['stimuli' filesep 'vowelA_stat_BF511_RAW_pos.wav']);
fsSig= 20e3;
sig_raw= helper.gen_resample(sig_raw, fsOrg, fsSig);
helper.plot_snap_fft_stat_vowel(uRatePos, uRateNeg, anl.fs, sig_raw, fsSig, anl.tStart, anl.tEnd, cur_unit_data.TC, 1);

if saveFig
    fName_psd= [DirStruct.latexDir 'Fig8'];
    saveas(gcf, fName_psd, 'epsc');
end