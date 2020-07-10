% function Fig2_corrGram_review(saveFig, LatexDir)
function Fig2_corrGram_review(saveFig, LatexDir)

if ~exist('saveFig', 'var')
    saveFig= 0;
end
if ~exist('LatexDir', 'var')
    LatexDir= ['figures' filesep];
end

figSize_cm= [15 5 19.05 8.5];
figure_prop_name = {'PaperPositionMode', 'units', 'Position', 'Renderer'};
figure_prop_val =  { 'auto', 'centimeters', figSize_cm, 'painters'};  % [Xcorner Ycorner Xwidth Ywidth]
figHan= 1;
figure(figHan);
clf;
set(gcf,figure_prop_name, figure_prop_val);

codesDir= pwd;
addpath(codesDir);

dataDir= ['data' filesep 'SK-2007_12_11-AN_normal' filesep];
fName= 'p0033_u1_07_SAM.m';
cd(dataDir);
cur_data= helper.loadpic(helper.getPicNum(fName));
target_dBSPL= 50;

all_calib_files= dir('*calib*');
calib_data= helper.loadpic(helper.getPicNum(all_calib_files(end).name));
calib_data= calib_data.CalibData;

pic_track_unit_NUMS= sscanf(fName, 'p%d_u%d_%d_*');
run(sprintf('Unit_%d_%02d', pic_track_unit_NUMS(2:3)));
unit_data= ans;
thresh= unit_data.Th;

fs= 20e3;
CF_Hz= 1e3*cur_data.Stimuli.condition.Carrfreq;
fm_Hz= 1e3*cur_data.Stimuli.condition.Modfrequency;
nCycles_fm= 2;
stimON= cur_data.Hardware.Trigger.StmOn/1e3;

Period_hist_dur= nCycles_fm/fm_Hz;
binEdges= 0:1/fs:Period_hist_dur;
binCenters= 1e3*(binEdges(1:end-1) + binEdges(2:end))/2;

calibOutatCF= interp1(calib_data(:,1), calib_data(:,2), CF_Hz/1e3);

% Get RLF and then SR
rlf_fName= dir(sprintf('p*u%d_%02d_RLV*', pic_track_unit_NUMS(2:3)));
rlf_data= helper.loadpic(helper.getPicNum(rlf_fName.name)); % p0021_u1_06_RLV

cd(codesDir);
rmpath(codesDir);

rlf_spls= calibOutatCF - rlf_data.Line.attens.Tone(:,2);
rlf_binEdges= (min(rlf_data.spikes{1}(:,1))-0.5):(max(rlf_data.spikes{1}(:,1))+0.5);
rlf_val= histcounts(rlf_data.spikes{1}(:,1), rlf_binEdges);
RLVparams= helper.fitRLfun(rlf_spls(:), rlf_val(:), 0, 0);
fprintf('SR for the SAM-example unit = %.1f\n', RLVparams.R_SP);


dBSPLs= calibOutatCF - cur_data.Line.attens.list(:,2);
unique_dBSPLs= unique(dBSPLs);
cur_dBSPL= unique_dBSPLs(abs(unique_dBSPLs-target_dBSPL)<2.6);
valid_lines= find(ismember(dBSPLs, cur_dBSPL));
curSpikes= cur_data.spikes{1};
abs_ref_period= 0.6e-3;

spike_trains= cell(length(valid_lines), 1);
for trainVar= 1:length(valid_lines)
    valid_spikes= ismember(curSpikes(:,1), valid_lines(trainVar));
    temp_spikes= curSpikes(valid_spikes, 2);
    temp_spikes= temp_spikes(temp_spikes<stimON);
    diff_isi= diff(temp_spikes);
    violation_index= 1 + find(diff_isi<abs_ref_period);
    temp_spikes(violation_index)= [];
    spike_trains{trainVar}= temp_spikes;
end

%% Compute correlograms
[sig, tSig]= helper.create_SAM(CF_Hz, fm_Hz, fs, 1, stimON, [], [], -pi/3);
sig= sig .* (sig>0);
sig= .1*sig/max(sig);
acf_sig= xcorr(sig);
delay_sig= (-(length(sig)-1):(length(sig)-1))/fs*1e3;


delayBinWidth= 1/fs;

[foISIhist,delay_fo] = helper.first_order_autoCorrGram(spike_trains, delayBinWidth, stimON);
[aoISIhist, delay_ao]= helper.all_order_autoCorrGram(spike_trains, delayBinWidth, stimON);
[shuf_auto_corrGram, delay_sac,AVGrate]= helper.SACfull_m(spike_trains, delayBinWidth, stimON);
NUMspikeREPS= length(spike_trains);
shuf_auto_corrGram= shuf_auto_corrGram*(NUMspikeREPS*(NUMspikeREPS-1)*AVGrate^2*delayBinWidth*stimON); % don't need normalization


%% Plot
yRange= 60;
nSProws= 2;
nSPcols= 4;
plt.xtick_CF_val= [50 1e3 4e3];
plt.xtick_CF_lab= cellfun(@(x) num2str(x), num2cell(plt.xtick_CF_val/1e3), 'uniformoutput', false);


% Time ACF
ax(1)= subplot(nSProws, nSPcols, 1);
ACFhan(1)= plot(delay_sig, acf_sig);
ttlHan(1)= title('Stimulus');
% xlabel('Delay (\mus)');
ylHan(1)= ylabel('Number of intervals', 'Units', 'normalized');

ax(2)= subplot(nSProws, nSPcols, 2);
ACFhan(2)= plot(delay_fo*1e3, foISIhist);
ttlHan(2)= title('FO-ISIH');
xlab_t_han= xlabel('Delay (ms)');

ax(3)= subplot(nSProws, nSPcols, 3);
ACFhan(3)= plot(delay_ao*1e3, aoISIhist);
ttlHan(3)= title('AO-ISIH');
% xlabel('Delay (\mus)');

ax(4)= subplot(nSProws, nSPcols, 4);
ACFhan(4)= plot(delay_sac*1e3, shuf_auto_corrGram);
ttlHan(4)= title('SAC');
% xlabel('Delay (\mus)');

linkaxes(ax, 'x');
xlim(ax(1), [-1 26]);
xlab_t_han.Position(1)= 30;

% ACF PSD
bx(1)= subplot(nSProws, nSPcols, 5);
[Pxx, ~, DFThan(1)]= helper.plot_dft(acf_sig, fs, 'yrange', yRange);
hold on;
text(fm_Hz, max(Pxx), 'F_m', 'HorizontalAlignment', 'center');
text(CF_Hz, max(Pxx)+7, 'F_c', 'HorizontalAlignment', 'center');
text(2*CF_Hz, max(Pxx)-10, '2F_c', 'HorizontalAlignment', 'center');
xlabel('');
ylHan(2)= ylabel('DFT magnitude (dB)', 'Units', 'normalized');


bx(2)= subplot(nSProws, nSPcols, 6);
[Pxx, ~, DFThan(2)]= helper.plot_dft(foISIhist, fs, 'yrange', yRange);
hold on;
text(fm_Hz, max(Pxx)-5, 'F_m', 'HorizontalAlignment', 'center');
text(CF_Hz, max(Pxx)-4, 'F_c', 'HorizontalAlignment', 'center');
ylabel('');
xlab_han= xlabel('Frequency (kHz)');

bx(3)= subplot(nSProws, nSPcols, 7);
[Pxx, ~, DFThan(3)]= helper.plot_dft(aoISIhist, fs, 'yrange', yRange);
hold on;
text(fm_Hz, max(Pxx), 'F_m', 'HorizontalAlignment', 'center');
text(2.5*fm_Hz, max(Pxx)-15, '2F_m', 'HorizontalAlignment', 'center');
text(CF_Hz, max(Pxx)+7, 'F_c', 'HorizontalAlignment', 'center');
text(2*CF_Hz, max(Pxx)-10, '2F_c', 'HorizontalAlignment', 'center');
ylabel('');
xlabel('');

bx(4)= subplot(nSProws, nSPcols, 8);
[Pxx, ~, DFThan(4)]= helper.plot_dft(shuf_auto_corrGram, fs, 'yrange', yRange);
hold on;
text(fm_Hz, max(Pxx), 'F_m', 'HorizontalAlignment', 'center');
text(2.5*fm_Hz, max(Pxx)-15, '2F_m', 'HorizontalAlignment', 'center');
text(CF_Hz, max(Pxx)+7, 'F_c', 'HorizontalAlignment', 'center');
text(2*CF_Hz, max(Pxx)-10, '2F_c', 'HorizontalAlignment', 'center');
ylabel('');
xlabel('');

set(bx, 'xtick', plt.xtick_CF_val, 'xticklabel', plt.xtick_CF_lab);
set(ACFhan, 'linew', 1, 'color', helper.get_color('b')); 
set(DFThan, 'linew', 1.5, 'color', helper.get_color('b')); 

linkaxes(bx, 'x');
xlim(bx(1), [fm_Hz/2 3*CF_Hz]);
xlab_han.Position(1)= 12e3;

set(findall(figHan,'-property','FontSize'),'FontSize', 9);
set(findall(figHan,'-property','Box'),'Box', 'off');
txtHan= helper.add_subplot_letter_txp(nSProws, nSPcols, 11, .05, 1.05);
set(ttlHan, 'FontSize', 11);

ylHan(1).Position(1)= max([ylHan(1).Position(1), ylHan(2).Position(1)]);
ylHan(2).Position(1)= max([ylHan(1).Position(1), ylHan(2).Position(1)]);

%% Set axes sizes
Xcorner_C= .06;
Xwidth_C= .2025;
Xshift_horz= .04;
Ycorner_C= .10;
Ywidth_C= .33;
Yshift_C= .17;

% A
set(ax(1),'Position',[Xcorner_C, Ycorner_C+Ywidth_C+Yshift_C, Xwidth_C, Ywidth_C])
drawnow

% B
set(ax(2),'Position',[Xcorner_C+Xwidth_C+Xshift_horz, Ycorner_C+Ywidth_C+Yshift_C, Xwidth_C, Ywidth_C])
drawnow

% C
set(ax(3),'Position',[Xcorner_C+2*Xwidth_C+2*Xshift_horz, Ycorner_C+Ywidth_C+Yshift_C, Xwidth_C, Ywidth_C])
drawnow

% D
set(ax(4),'Position',[Xcorner_C+3*Xwidth_C+3*Xshift_horz, Ycorner_C+Ywidth_C+Yshift_C, Xwidth_C, Ywidth_C])
drawnow

% E
set(bx(1),'Position',[Xcorner_C, Ycorner_C, Xwidth_C, Ywidth_C])
drawnow

% F
set(bx(2),'Position',[Xcorner_C+Xwidth_C+Xshift_horz, Ycorner_C, Xwidth_C, Ywidth_C])
drawnow

% G
set(bx(3),'Position',[Xcorner_C+2*Xwidth_C+2*Xshift_horz, Ycorner_C, Xwidth_C, Ywidth_C])
drawnow

% H
set(bx(4),'Position',[Xcorner_C+3*Xwidth_C+3*Xshift_horz, Ycorner_C, Xwidth_C, Ywidth_C])
drawnow
%%

if saveFig
    saveas(gcf, [LatexDir 'Fig2'], 'epsc');
end
