% function Fig1_intro_nonstat_demo(saveFig, LatexDir)
function Fig1_intro_nonstat_demo(saveFig, LatexDir)

if ~exist('saveFig', 'var')
    saveFig= 0;
end
if ~exist('LatexDir', 'var')
    LatexDir= ['figures' filesep];
end

figHan= 1;


%% Load speech
anl.speech.tStart= 0.25; %1.55;
anl.speech.tEnd= 2.35; %2.35;
anl.speech.fs= 12e3;

CodesDir= pwd;

[sig_eng, fs_eng]= audioread(['stimuli' filesep 'Stim_S_P.wav']);
sig_eng= helper.gen_resample(sig_eng, fs_eng, anl.speech.fs);
fs_eng= anl.speech.fs;
t_eng= (1:length(sig_eng))/fs_eng;
anl.speech.tMask= t_eng>anl.speech.tStart & t_eng<anl.speech.tEnd;

[pSPL, tSPL]= helper.gen_get_spl_vals(sig_eng, fs_eng, 20e-3, .8);
validInds= tSPL>anl.speech.tStart & tSPL<anl.speech.tEnd;

%
eng_data.voiced_boundaries= helper.find_voicing_boundaries(sig_eng, fs_eng, 0, .13);
eng_data.voiced_boundaries(:,1)= eng_data.voiced_boundaries(:,1)+10e-3;
eng_data.voiced_boundaries(:,2)= eng_data.voiced_boundaries(:,2)-10e-3;
eng_data.voiced_inds= any(t_eng(anl.speech.tMask)>eng_data.voiced_boundaries(:,1) & t_eng(anl.speech.tMask)<eng_data.voiced_boundaries(:,2), 1);

temp_formants= load(['data' filesep 'english_formant.mat']);
temp_formants= temp_formants.formant_data;
eng_data.trajectory.f1= interp1([temp_formants.time], [temp_formants.f1], t_eng(anl.speech.tMask), 'pchip')/1e3;
eng_data.trajectory.f1(~eng_data.voiced_inds)= nan;
eng_data.trajectory.f2= interp1([temp_formants.time], [temp_formants.f2], t_eng(anl.speech.tMask), 'pchip')/1e3;
eng_data.trajectory.f2(~eng_data.voiced_inds)= nan;
eng_data.trajectory.f3= interp1([temp_formants.time], [temp_formants.f3], t_eng(anl.speech.tMask), 'pchip')/1e3;
eng_data.trajectory.f3(~eng_data.voiced_inds)= nan;

% Start plotting
plt.lw3= 1.5;
plt.nSProws= 4;
plt.nSPcols= 3;
anl.speech.tWindow= 60e-3;
plt.spect.fracOVlap= .9;
plt.spect.nfft= 2^(1+nextpow2(anl.speech.fs*anl.speech.tWindow));
plt.spect.useDefaultPlot= false;
plt.spect.tStart= anl.speech.tStart + anl.speech.tWindow/2;

%% 
figSize_cm= [15 5 19.05 10];
figure_prop_name = {'PaperPositionMode','units','Position', 'Renderer'};
figure_prop_val =  { 'auto'            ,'centimeters', figSize_cm, 'painters'};  % [Xcorner Ycorner Xwidth Ywidth]
figure(figHan);
clf;
set(gcf,figure_prop_name,figure_prop_val);

ax(1)= subplot(plt.nSProws, plt.nSPcols, [1 4 7]);
ax(2)= subplot(plt.nSProws, plt.nSPcols, 10);
ax(3)= subplot(plt.nSProws, plt.nSPcols, [2 5]);
ax(4)= subplot(plt.nSProws, plt.nSPcols, [8 11]);
ax(5)= subplot(plt.nSProws, plt.nSPcols, [3 6]);
ax(6)= subplot(plt.nSProws, plt.nSPcols, [9 12]);

%
Xcorner_X= .05;
Xwidth_X= .26;
Xshift_X= .08;

Ycorner_X= .09;
Ywidth_A= .48;
Ywidth_B= .28;
Ywidth_X= .38;
Yshift_X= .1;

% B
set(ax(2),'Position',[Xcorner_X, Ycorner_X, Xwidth_X, Ywidth_B])
drawnow

        % A
set(ax(1),'Position',[Xcorner_X, Ycorner_X+Ywidth_B+Yshift_X, Xwidth_X, Ywidth_A])
drawnow

% D
set(ax(4),'Position',[Xcorner_X+Xwidth_X+1.5*Xshift_X, Ycorner_X, Xwidth_X, Ywidth_X])
drawnow

% C
set(ax(3),'Position',[Xcorner_X+Xwidth_X+1.5*Xshift_X, Ycorner_X+Ywidth_X+Yshift_X, Xwidth_X, Ywidth_X])
drawnow

% F
set(ax(6),'Position',[Xcorner_X+2*Xwidth_X+2*Xshift_X, Ycorner_X, Xwidth_X, Ywidth_X])
drawnow

% E
set(ax(5),'Position',[Xcorner_X+2*Xwidth_X+2*Xshift_X, Ycorner_X+Ywidth_X+Yshift_X, Xwidth_X, Ywidth_X])
drawnow

%%
axes(ax(2)); 
plot(t_eng(anl.speech.tMask)*1e3, sig_eng(anl.speech.tMask), 'color', helper.get_color('b'));
ylabHan(1)= ylabel('Amplitude', 'Units', 'normalized');
txt(2)= text(.05, 1.05, 'B', 'Units', 'normalized');

yyaxis right;
plot(tSPL(validInds)*1e3, pSPL(validInds), 'LineWidth', plt.lw3);
box off;
ylabel('Intensity (dB SPL)');

axes(ax(1));
hold on;
% timePlot= helper.plot_spectrogram(sig_eng(anl.speech.tMask), fs_eng, anl.speech.tWindow, plt.spect.fracOVlap, plt.spect.nfft, plt.spect.useDefaultPlot, plt.spect.tStart);
plt.ColBar= max(get(colorbar, 'Limits'));
caxis([plt.ColBar-40 plt.ColBar]);
colorbar off;
ylabHan(2)= ylabel('Frequency (kHz)', 'Units', 'normalized');
tFormant= 1e3*(anl.speech.tWindow/2+t_eng(anl.speech.tMask));
% tFormant(tFormant<min(timePlot))= nan;

lformant_Han(1)= line(tFormant, eng_data.trajectory.f1, plt.ColBar*ones(size(eng_data.trajectory.f1)),'linew', plt.lw3);
lformant_Han(2)= line(tFormant, eng_data.trajectory.f2, plt.ColBar*ones(size(eng_data.trajectory.f1)), 'linew', plt.lw3);
lformant_Han(3)= line(tFormant, eng_data.trajectory.f3, plt.ColBar*ones(size(eng_data.trajectory.f1)), 'linew', plt.lw3);
set(lformant_Han, 'color', 'k');
txt(1)= text(.05, 1.05, 'A', 'Units', 'normalized');

linkaxes(ax(1:2), 'x');
xlim([anl.speech.tStart anl.speech.tEnd]*1e3);
ylabHan(1).Position(1)= mean([ylabHan(1).Position(1), ylabHan(2).Position(1)]);
ylabHan(2).Position(1)= mean([ylabHan(1).Position(1), ylabHan(2).Position(1)]);

%% Load PSTH to tone
plt.lw= 1.5;

tone_files.tonePST= ['data' filesep 'SP-2016_08_08-Q245_AN' filesep 'p0029_u1_05_PST.mat'];
tone_files.SR= ['data' filesep 'SP-2016_08_08-Q245_AN' filesep 'p0028_u1_05_SR.mat'];

temp_data= load(tone_files.tonePST);
temp_data= temp_data .data;


sr_data= load(tone_files.SR);
sr_data= sr_data.data;
sr_spikes= sr_data.spikes{1};
sr_value= size(sr_spikes,1)/max(sr_spikes(:,1));
fprintf('SR for the tone-example unit = %.1f\n', sr_value);

fs= 20e3;
stim_dur= temp_data.Hardware.Trigger.StmOn/1e3;
binEdges_perHist= 0:1/fs:stim_dur;
fc= temp_data.Stimuli.main.tone.freq;
temp_spike_data= temp_data.spikes{1};
temp_spike_data= temp_spike_data(:,2);

xx= load(['data' filesep 'SP-2016_08_08-Q245_AN' filesep 'p0006_calib.mat']);
db_spl_used= interp1(xx.data.CalibData(:,1), xx.data.CalibData(:,2), fc/1e3) - unique(temp_data.Line.attens.Tone(:,2));
fprintf('Intensity for tone = %.1f dB SPL \n', db_spl_used);


t_uRate= (binEdges_perHist(1:end-1)+binEdges_perHist(2:end))/2;
uRate_pos= histcounts(temp_spike_data, binEdges_perHist);

fsPerHist= 40e3;
tStart= 20e-3;
binEdges_PerHist= 0:1/fsPerHist:1/fc;
perHistSpikeTimes= rem(temp_spike_data(temp_spike_data>tStart), 1/fc); 
PerHist_pos= histcounts(perHistSpikeTimes, binEdges_PerHist);
PerHist_pos= circshift(PerHist_pos, -15); % Just to make it look aligned with the +ve half-cycle
t_PerHist= (binEdges_PerHist(1:end-1)+binEdges_PerHist(2:end))/2;

perHist_fs= 100e3; 
[perHist_sig, perHist_time]= helper.create_sinusoid(fc, perHist_fs, 1/fc, 1, 0);

env_params.tStart= 7e-3;
env_params.tEnd= 100e-3;
[peaks, locs]= findpeaks(uRate_pos, 'MinPeakDistance', round(fs/2/fc),'MinPeakHeight', 5);
valid_loc_inds= (t_uRate(locs)>env_params.tStart) & (t_uRate(locs)<env_params.tEnd);
locs= locs(valid_loc_inds);

tone_stim= helper.gen_ramp_signal(sin(2*pi*fc*t_uRate), fs, 5e-3);
shift_tone= 6e-3;
tone_stim= [zeros(1, round(shift_tone*fs)), tone_stim];

anl.tone.tStart= 0e-3;
anl.tone.tEnd= 46.5e-3;
anl.tone.tMask= t_uRate>anl.tone.tStart & t_uRate<anl.tone.tEnd;

axes(ax(3));
hold on;
plot(t_uRate(anl.tone.tMask)*1e3, uRate_pos(anl.tone.tMask), 'color', helper.get_color('gray'), 'linew', plt.lw);
plot(t_uRate(anl.tone.tMask)*1e3, -2.5+2*tone_stim(anl.tone.tMask), 'color', helper.get_color('b'), 'linew', plt.lw);

xlim([anl.tone.tStart anl.tone.tEnd]*1e3);
txt(3)= text(.05, 1.05, 'C', 'Units', 'normalized');
ylim([-5 25]);


axes(ax(4));
hold on;
plot(t_PerHist*1e3, PerHist_pos, 'color', helper.get_color('gray'), 'linew', plt.lw);
plot(perHist_time*1e3, min(PerHist_pos) -50 + 75*perHist_sig, 'color', helper.get_color('b'), 'linew', plt.lw);
% ylabel('Number of spikes');
txt(5)= text(.05, 1.05, 'D', 'Units', 'normalized');
% legend('Per. Hist.', 'Stimulus', 'box', 'off', 'Location', 'northeast');
xlab_han= xlabel('Time (ms)');
ylabHan= ylabel('Discharge rate (spikes/bin)');
ylabHan.Position(1)= -0.27;
ylabHan.Position(2)= 1000;


%% Add codesDir to path
addpath(CodesDir);

sam_files.dataDir= ['data' filesep 'SK-2007_08_23-AN_normal' filesep];
cd(sam_files.dataDir);

sam_files.fName= 'p0029_u1_07_SAM.m';
cur_data= helper.loadpic(helper.getPicNum(sam_files.fName));
target_dBSPL= 35;

sam_files.calib_filename= 'p0003_calib.m';
calib_data= helper.loadpic(helper.getPicNum(sam_files.calib_filename));
calib_data= calib_data.CalibData;

sam_files.unit_filename= 'Unit_1_07';
run(sam_files.unit_filename);

% Get RLF and then SR
sam_files.rlf_filename= 'p0027_u1_07_RLV';
rlf_data= helper.loadpic(helper.getPicNum(sam_files.rlf_filename)); % p0021_u1_06_RLV

%% Remove codesDir from path
unit_data= ans;
thresh= unit_data.Th;

CF_Hz= 1e3*cur_data.Stimuli.condition.Carrfreq;
fm_Hz= 1e3*cur_data.Stimuli.condition.Modfrequency;
nCycles_fm= 1;
stimON= cur_data.Hardware.Trigger.StmOn/1e3;


tStart_psth= 237.5e-3;
tEnd_psth= 315e-3;
binEdges_PSTH= tStart_psth:1/fs:tEnd_psth;
binCenters_PSTH= 1e3*(binEdges_PSTH(1:end-1) + binEdges_PSTH(2:end))/2;
binCenters_PSTH= binCenters_PSTH-min(binCenters_PSTH);


Period_hist_dur= nCycles_fm/fm_Hz;
binEdges_perHist= 0:1/fs:Period_hist_dur;
binCenters_perHist= 1e3*(binEdges_perHist(1:end-1) + binEdges_perHist(2:end))/2;

calibOutatCF= interp1(calib_data(:,1), calib_data(:,2), CF_Hz/1e3);


rlf_spls= calibOutatCF - rlf_data.Line.attens.Tone(:,2);
rlf_binEdges= (min(rlf_data.spikes{1}(:,1))-0.5):(max(rlf_data.spikes{1}(:,1))+0.5);
rlf_val= histcounts(rlf_data.spikes{1}(:,1), rlf_binEdges);
RLVparams= helper.fitRLfun(rlf_spls(:), rlf_val(:), 0, 0);
fprintf('SR for the SAM-example unit = %.1f\n', RLVparams.R_SP);


dBSPLs= calibOutatCF - cur_data.Line.attens.list(:,2);
unique_dBSPLs= unique(dBSPLs);
cur_dBSPL= unique_dBSPLs(abs(unique_dBSPLs-target_dBSPL)<5);
valid_lines= find(ismember(dBSPLs, cur_dBSPL));
curSpikes= cur_data.spikes{1};
valid_spikes= ismember(curSpikes(:,1), valid_lines);
spike_times= curSpikes(valid_spikes, 2);
spike_times= spike_times(spike_times<stimON);


uRate_PSTH= histcounts(spike_times, binEdges_PSTH);

spike_times_perHist= rem(spike_times, Period_hist_dur);
uRate_perHist= histcounts(spike_times_perHist, binEdges_perHist);
[sig_perHist, tSig_perHist]= helper.create_SAM(CF_Hz, fm_Hz, fs, 1, nCycles_fm/fm_Hz, [], [], -pi/3);
samAmp_perHist= 3.5;

[sig_PSTH, tSig_PSTH]= helper.create_SAM(CF_Hz, fm_Hz, fs, 1, tEnd_psth-tStart_psth, [], [], -pi/3);
samAmp_PSTH= 0.75;

cd(CodesDir);
rmpath(CodesDir);

axes(ax(5));
hold on;
plot(binCenters_PSTH, uRate_PSTH, 'color', helper.get_color('gray'), 'linew', plt.lw);
plot(1e3*tSig_PSTH, -samAmp_PSTH+ .75*samAmp_PSTH*.75*sig_PSTH/max(sig_PSTH), 'color', helper.get_color('b'), 'linew', plt.lw);
plot(1e3*tSig_PSTH, -samAmp_PSTH+ abs(hilbert(.75*samAmp_PSTH*.75*sig_PSTH/max(sig_PSTH))), 'color', helper.get_color('r'), 'linew', plt.lw);
txt(4)= text(.05, 1.05, 'E', 'Units', 'normalized');
ylim([-1.3 2.1]);

% plot
axes(ax(6));
hold on;
plot(binCenters_perHist, uRate_perHist, 'color', helper.get_color('gray'), 'linew', plt.lw);
plot(1e3*tSig_perHist, -1.2*samAmp_perHist + samAmp_perHist*sig_perHist/max(sig_perHist), 'color', helper.get_color('b'), 'linew', plt.lw);
plot(1e3*tSig_perHist, -1.2*samAmp_perHist + abs(hilbert(samAmp_perHist*sig_perHist/max(sig_perHist))), 'color', helper.get_color('r'), 'linew', plt.lw);

% ylabel('Number of spikes');
% xlabel('Time (ms)');

txt(6)= text(.05, 1.05, 'F', 'Units', 'normalized');
xlim([0 1e3*(nCycles_fm + 0.05)/fm_Hz]);
ylim([-2.3*samAmp_perHist max(uRate_perHist)+0.5]);

set(findall(figHan,'-property','FontSize'),'FontSize', 9);
set(txt, 'FontSize', 11);


if saveFig
    print([LatexDir 'Fig1'], '-dpng',  '-r600');
    print([LatexDir 'Fig1'], '-dtiff',  '-r600');
end