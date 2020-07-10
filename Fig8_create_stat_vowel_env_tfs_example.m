% This script estimates PSDs for individual units (for both statioanry and
% kinematic condition) and saves it in a folder,

clear;
clc;

saveFig= 0;
% Init params
figHan.time= 1;
figHan.psd= 2;

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
DirStruct.INdata= ['data' filesep];

ChinID= 374;
SynCapData= load(sprintf('%sQ%d_allSyncCap.mat', DirStruct.INdata, ChinID));
cur_unit_data=SynCapData.SynCapData(4);


uRatePos= histcounts(cell2mat(cur_unit_data.stat.RAW.pos'), anl.tEdge_hist);
uRateNeg= histcounts(cell2mat(cur_unit_data.stat.RAW.neg'), anl.tEdge_hist);

[sig_raw, fsOrg]= audioread(['stimuli' filesep 'vowelA_stat_BF511_RAW_pos.wav']);
fsSig= 20e3;
sig_raw= helper.gen_resample(sig_raw, fsOrg, fsSig);
PSDout= helper.plot_snap_fft_stat_vowel(uRatePos, uRateNeg, anl.fs, sig_raw, fsSig, anl.tStart, anl.tEnd, cur_unit_data.TC, 1);

if saveFig
    DirStruct.OUT_fig_psd_single= '/home/parida/Dropbox/Articles/neural_temporal_coding/figures/';
    fName_psd= [DirStruct.OUT_fig_psd_single 'Fig8'];
    saveas(gcf, fName_psd, 'epsc');
end