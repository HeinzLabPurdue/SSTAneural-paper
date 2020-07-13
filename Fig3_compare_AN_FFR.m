% function Fig3_compare_AN_FFR(saveFig, LatexDir)
function Fig3_compare_AN_FFR(saveFig, LatexDir)

if ~exist('saveFig', 'var')
    saveFig= 0;
end
if ~exist('LatexDir', 'var')
    LatexDir= ['figures' filesep];
end

figSize_cm= [15 5 13.2 8];
figure_prop_name = {'PaperPositionMode','units','Position', 'Renderer'};
figure_prop_val =  { 'auto'            ,'centimeters', figSize_cm, 'painters'};  % [Xcorner Ycorner Xwidth Ywidth]
figure(1);
clf;
set(gcf,figure_prop_name,figure_prop_val);

ffr.dataDir= ['data' filesep 'SP-2019_06_09-Q374_SFR_pink500Hz_NH' filesep];
% For 374, two files for clean speech. So should use the combined
% file. That's why just 3 files.

ANdataDir= ['data' filesep];



%% First load FFR data
ffr.file= dir([ffr.dataDir 'a*SSN_Stim_S_P*']);
ffr.struct_data= load([ffr.file.folder filesep ffr.file.name]);
ffr.struct_data= ffr.struct_data.data;

ffr.fs= ffr.struct_data.Stimuli.RPsamprate_Hz;
ffr.data.pos= ffr.struct_data.AD_Data.AD_Avg_PO_V{1};
ffr.data.neg= ffr.struct_data.AD_Data.AD_Avg_NP_V{1};
ffr.data.tfs= (ffr.data.pos-ffr.data.neg)/2;
ffr.data.env= (ffr.data.pos+ffr.data.neg)/2;

%% Load all AN data for Q374
chinDanishData= [];

Chin_Q374= 374;
an.all_files= dir([ANdataDir '*DirBased*']);

SynCapData= load(sprintf('%sQ%d_allSyncCap.mat', ANdataDir, Chin_Q374));
SynCapData=SynCapData.SynCapData;

for unitVar= 1:length(SynCapData)
    cur_unit_data= SynCapData(unitVar);
    
    cur_delay= 4.59e-3; %This is the additional delay of using RP2-inv-calib

    chinDanishData(unitVar).CF_Hz= cur_unit_data.CF_Hz;  %#ok<*AGROW>
    chinDanishData(unitVar).SR= cur_unit_data.SR;
    chinDanishData(unitVar).pos= cell2mat(cur_unit_data.danish.pos)-cur_delay;
    chinDanishData(unitVar).neg= cell2mat(cur_unit_data.danish.neg)-cur_delay;
end

%% All other NH chins
dirStruct.loading_dir= ['data' filesep 'DanishData' filesep];
all_chinID= [321 322 325 338 341 343 346 347 354 355 373 379];

allChinSpikeData = [];
if ~exist('all_chinID', 'var')
    allfiles=dir([dirStruct.loading_dir '*.mat']);
else
    for fileVar=1:length(all_chinID)
        curChinID= all_chinID(fileVar);
        allfiles(fileVar)=dir([dirStruct.loading_dir '*' num2str(curChinID) '*']);
    end
end
for fileVar=1:length(allfiles)
    temp = load([dirStruct.loading_dir allfiles(fileVar).name]);
    allChinSpikeData = [allChinSpikeData; temp.spike_data']; 
end

[~, uniqInds]= unique([ [allChinSpikeData.chinID]', [allChinSpikeData.track]', [allChinSpikeData.unit]', [allChinSpikeData.SPL]'], 'rows');
uniqSpikeData= allChinSpikeData(uniqInds);
clear allChinSpikeData;

%% Combine Q374 and other chins
count= length(chinDanishData);

for unitVar= 1:length(uniqSpikeData)
    count= count+1;
    
    
    if uniqSpikeData(unitVar).chinID>368 % For newer animals
        cur_delay= 4.59e-3;
    else
        cur_delay= 0;
    end

    chinDanishData(count).CF_Hz= uniqSpikeData(unitVar).CF_Hz;
    chinDanishData(count).SR= uniqSpikeData(unitVar).SR;
    chinDanishData(count).pos= cell2mat(uniqSpikeData(unitVar).SpikeTrains{1,1})-cur_delay;
    chinDanishData(count).neg= cell2mat(uniqSpikeData(unitVar).SpikeTrains{1,2})-cur_delay;
end

%% Construct PSTH
anl.stimPeriod= 1800e-3; % stim rep rate (stim on + stim off)
anl.stimDuration= 1300e-3;
anl.fs= 10e3; % 1/(bin resolution) for uRate/histcount
anl.nw= 1.5;
anl.tEdge_hist= 0:1/anl.fs:anl.stimDuration;
anl.tBinCenter= (anl.tEdge_hist(1:end-1) + anl.tEdge_hist(2:end))/2;

all_an.psth_compound= nan(length(chinDanishData), length(anl.tBinCenter));
all_an.psth_sum= nan(length(chinDanishData), length(anl.tBinCenter));
all_an.psth_pos= nan(length(chinDanishData), length(anl.tBinCenter));
all_an.psth_neg= nan(length(chinDanishData), length(anl.tBinCenter));

for unitVar= 1:length(chinDanishData)
    cur_unit_data= chinDanishData(unitVar);
    if ~isempty(cur_unit_data.pos)
        temp_danish_uRate_pos= histcounts(cur_unit_data.pos, anl.tEdge_hist)/length(cur_unit_data.pos);
        temp_danish_uRate_neg= histcounts(cur_unit_data.neg, anl.tEdge_hist)/length(cur_unit_data.neg);
        all_an.psth_pos(unitVar, :)= temp_danish_uRate_pos;
        all_an.psth_neg(unitVar, :)= temp_danish_uRate_neg;
        all_an.psth_compound(unitVar, :)= (temp_danish_uRate_pos - temp_danish_uRate_neg)/2;
        all_an.psth_sum(unitVar, :)= (temp_danish_uRate_pos + temp_danish_uRate_neg)/2;
    end
end

%% Apply equal weight to CF normalize for CF distribution
allCFs= [chinDanishData.CF_Hz];
CF_edges= logspace(log10(196), log10(11e3), 18); % ~one-third ocave
CF_bins= ( CF_edges(1:end-1) + CF_edges(2:end)) / 2;

an.psth_compound= nan(length(CF_bins), length(anl.tBinCenter));
an.psth_sum= nan(length(CF_bins), length(anl.tBinCenter));
an.psth_pos= nan(length(CF_bins), length(anl.tBinCenter));
an.psth_neg= nan(length(CF_bins), length(anl.tBinCenter));

for binVar= 1:length(CF_bins)
    CFlow= CF_edges(binVar);
    CFhi= CF_edges(binVar+1);
    
    validCFs= (allCFs > CFlow) & (allCFs < CFhi);
    an.psth_pos(binVar, :)= nanmean(all_an.psth_pos(validCFs, :));
    an.psth_neg(binVar, :)= nanmean(all_an.psth_neg(validCFs, :));
    an.psth_compound(binVar, :)= nanmean(all_an.psth_compound(validCFs, :));
    an.psth_sum(binVar, :)= nanmean(all_an.psth_sum(validCFs, :));
end

%% Load stim 
fs= 10e3;

[sig, fs_sig]= audioread(['stimuli' filesep 'FLN_Stim_S_P.wav']);
sig= helper.gen_resample(sig, fs_sig, fs);
t_sig= (1:length(sig))/fs;


%% Plot 
figure(1);
clf;
plt.lw1= 0.62;
plt.lw2= 1.0;
plt.lw3= 1.25;
plt.tick_len= [.025 .025];
xtick_vals_freq= [100 500 1e3];
xtick_labs_freq= cellfun(@(x) num2str(x), num2cell(xtick_vals_freq/1e3), 'UniformOutput', false);

% FFR
Hd_BP= helper.get_filter_designfilt('bp', [100 1.5e3], fs, 2);

ffr_tfs_plot= helper.gen_resample(ffr.data.tfs, ffr.fs, fs);
ffr_tfs_plot= filtfilt(Hd_BP, ffr_tfs_plot);
ffr_env_plot= helper.gen_resample(ffr.data.env, ffr.fs, fs);
ffr_env_plot= filtfilt(Hd_BP, ffr_env_plot);
ffr_time_plot= (1:length(ffr_tfs_plot))/fs;
ffr_tMask= ffr_time_plot<max(t_sig);

% AN PSTH
an_tfs_plot= helper.gen_resample(nanmean(all_an.psth_compound, 1), anl.fs, fs);
an_tfs_plot= filtfilt(Hd_BP, an_tfs_plot);
an_env_plot= helper.gen_resample(nanmean(all_an.psth_sum, 1), anl.fs, fs);
an_env_plot= filtfilt(Hd_BP, an_env_plot);
an_time_plot= (1:length(an_tfs_plot))/fs;


equal_rms_an_ffr= 1; 
if equal_rms_an_ffr
    ffr_tfs_plot= ffr_tfs_plot/rms(ffr_tfs_plot);
    ffr_env_plot= -ffr_env_plot/rms(ffr_env_plot);
    an_tfs_plot= an_tfs_plot/rms(an_tfs_plot);
    an_env_plot= an_env_plot/rms(an_env_plot);
else 
    ffr_tfs_plot= ffr_tfs_plot/max(abs(ffr_tfs_plot));
    an_tfs_plot= .75*an_tfs_plot/max(abs(an_tfs_plot));
    ffr_env_plot= -ffr_env_plot/max(abs(ffr_env_plot));
    an_env_plot= an_env_plot/max(abs(an_env_plot));
end

%%
clf;

tStart= 0.23; tEnd= 0.33;  % Perfect
F0_freq= 100;
F1_freq= 670; 
F2_freq= 1410; 
% tStart= 1.05; tEnd= 1.2; % Strong F1, not super helpful

da_time= 20e-3;
ga_time= 200e-3;
s_time= 740e-3;

ffr_valid_inds= ffr_time_plot>tStart & ffr_time_plot<tEnd;
an_valid_inds= an_time_plot>tStart & an_time_plot<tEnd;
sig_valid_inds= t_sig>tStart & t_sig<tEnd; 

A_Shift_TFS= 2*max(abs(an_tfs_plot));
A_Shift_ENV= 1.5*max(abs(an_env_plot));

ax(1)= subplot(2,5,1:3);
hold on;
plot(ffr_time_plot(ffr_tMask), A_Shift_TFS+ffr_tfs_plot(ffr_tMask), 'color', helper.get_color('b'));
plot(an_time_plot, an_tfs_plot, 'color', helper.get_color('r'));
plot(t_sig, -A_Shift_ENV/2+A_Shift_ENV*sig, 'color', helper.get_color('k'));
plot([tStart tStart], [-2*A_Shift_TFS 2*A_Shift_TFS], '--', 'color', helper.get_color('prp'), 'linew', plt.lw2);
plot([tEnd tEnd], [-2*A_Shift_TFS 2*A_Shift_TFS], '--', 'color', helper.get_color('prp'), 'linew', plt.lw2);
txtHan(1)= text(.05, 1.05, 'A. Difference', 'units', 'normalized');
set(gca, 'XTickLabel', '');
ylim([-1.75*A_Shift_TFS 1.75*A_Shift_TFS]);

bx(1)= subplot(2,5,4:5);
hold on;
[Pxx_stim,~, TFS_PSD_han(3)]= helper.plot_dpss_psd(3*A_Shift_ENV*sig(sig_valid_inds), fs, 'nw', anl.nw);
[~,~, TFS_PSD_han(1)]= helper.plot_dpss_psd(ffr_tfs_plot(ffr_valid_inds), fs, 'nw', anl.nw);
[~,~, TFS_PSD_han(2)]= helper.plot_dpss_psd(an_tfs_plot(an_valid_inds), fs, 'nw', anl.nw);
text(F0_freq, max(Pxx_stim)+1, 'F_0', 'HorizontalAlignment', 'center'); 
text(F1_freq, max(Pxx_stim)+4, 'F_1', 'HorizontalAlignment', 'center'); 
text(F2_freq, max(Pxx_stim)-3, 'F_2', 'HorizontalAlignment', 'center'); 

txtHan(2)= text(.05, 1.05, 'B. Difference spectrum', 'units', 'normalized');
xlabel('');
set(TFS_PSD_han(1), 'linew', plt.lw3, 'color', helper.get_color('b'));
set(TFS_PSD_han(2), 'linew', plt.lw2, 'LineStyle', '-', 'color', helper.get_color('r'));
set(TFS_PSD_han(3), 'linew', plt.lw3, 'LineStyle', '-', 'Color', helper.get_color('k'));
ylabel('');
set(gca, 'XTick', xtick_vals_freq, 'XTickLabel', '');
% set(gca, 'XTickLabel', '');

ax(2)= subplot(2,5,6:8);
hold on;
lineHan(1)= plot(ffr_time_plot(ffr_tMask), A_Shift_ENV+ffr_env_plot(ffr_tMask), 'color', helper.get_color('b'));
lineHan(2)= plot(an_time_plot, an_env_plot, 'color', helper.get_color('r'));
lineHan(3)= plot(nan, nan, 'k');
plot([tStart tStart], [-1.5*A_Shift_TFS 2.5*A_Shift_TFS], '--', 'color', helper.get_color('prp'), 'linew', plt.lw2);
plot([tEnd tEnd], [-1.5*A_Shift_TFS 2.5*A_Shift_TFS], '--', 'color', helper.get_color('prp'), 'linew', plt.lw2);
text(da_time, -0.7*A_Shift_TFS, '/d/', 'HorizontalAlignment', 'center');
text(ga_time, -0.7*A_Shift_TFS, '/g/', 'HorizontalAlignment', 'center');
text(s_time, -0.7*A_Shift_TFS, '/s/', 'HorizontalAlignment', 'center');
txtHan(3)= text(.05, 1.05, 'C. Sum', 'units', 'normalized');
xlabel('Time (s)')
yLab_time= ylabel('Normalized Amplitude');
yLab_time.Position(1)= -.18;
yLab_time.Position(2)= 1.6*A_Shift_ENV;
ylim(ax(2), [-.85*A_Shift_ENV 1.5*A_Shift_ENV]);
[lgHan, icons]= legend(lineHan, 'FFR', 'PSTH', 'SIG');

bx(2)= subplot(2,5,9:10);
hold on;
% [~,~, ENV_PSD_han(3)]= plot_dpss_psd(2*A_Shift_ENV*sig, fs, 'nw', anl.nw);
[~,~, ENV_PSD_han(1)]= helper.plot_dpss_psd(ffr_env_plot(ffr_valid_inds), fs, 'nw', anl.nw);
[~,~, ENV_PSD_han(2)]= helper.plot_dpss_psd(an_env_plot(an_valid_inds), fs, 'nw', anl.nw);
txtHan(4)= text(.05, 1.05, 'D. Sum spectrum', 'units', 'normalized');
set(gca, 'XTick', xtick_vals_freq, 'XTickLabel', xtick_labs_freq);
text(F0_freq, max(Pxx_stim), 'F_0', 'HorizontalAlignment', 'center'); 
xlabel('Frequency (kHz)');
set(ENV_PSD_han(1), 'linew', plt.lw3, 'color', helper.get_color('b'));
set(ENV_PSD_han(2), 'linew', plt.lw2, 'LineStyle', '-', 'color', helper.get_color('r'));
% set(ENV_PSD_han(3), 'linew', plt.lw3, 'LineStyle', '-', 'Color', helper.get_color('k'));

linkaxes(ax, 'x')
xlim(ax(1), [-0.02 1.45]);
% ylim([-1.5 2.5])


linkaxes(bx);
subplot(2,5,9:10);
xlim([75 1.75e3]);
ylim([-50 0]);
yLab_psd= ylabel('PSD (dB/Hz)');
yLab_psd.Position(1)= 47;
yLab_psd.Position(2)= max(ylim);

% Update legend locations
lgHan.Box= 'off'; 
lgHan.FontSize= 9;
lgHan.Position(1:2)= [.4 .1];

icons(4).XData= mean(icons(4).XData) + [0.1 .3];
icons(4).LineWidth= 1;
icons(6).XData= mean(icons(6).XData) + [0.1 .3];
icons(6).LineWidth= 1;
icons(8).XData= mean(icons(8).XData) + [0.1 .3];
icons(8).LineWidth= 1;

%% define new axes for AB
Xshift_horz= .1;

Xwidth_AC= .45;
Ywidth_AC= .39;
Xcorner_AC= .08;
Yshift_AC= .06;
Ycorner_AC= .115;

% B
set(ax(2),'Position',[Xcorner_AC Ycorner_AC Xwidth_AC Ywidth_AC])
drawnow
% A
set(ax(1),'Position',[Xcorner_AC Ycorner_AC+Ywidth_AC+Yshift_AC Xwidth_AC Ywidth_AC])
drawnow
%% define new axes for AB
Xwidth_BD= .36;
Ywidth_BD= Ywidth_AC;
Xcorner_BD= Xcorner_AC+Xwidth_AC+Xshift_horz;
Yshift_BD= Yshift_AC;
Ycorner_BD= Ycorner_AC;
% Yshift=0.06;

% D
set(bx(2),'Position',[Xcorner_BD Ycorner_BD Xwidth_BD Ywidth_BD])
drawnow
% C
set(bx(1),'Position',[Xcorner_BD Ycorner_BD+Ywidth_BD+Yshift_BD Xwidth_BD Ywidth_BD])
drawnow

set(findall(gcf,'-property','FontSize'),'FontSize', 9);
set(txtHan,'FontSize', 10);
set(findall(gcf,'-property','TickLength'),'TickLength', plt.tick_len, 'units', 'normalized');


% fName= 'FFR_AN_compare';
fName= 'Fig3';
if saveFig
    saveas(gcf, [LatexDir fName], 'epsc');
end
fprintf('total fibers (N)=%d\n', sum(cellfun(@(x) ~isempty(x), {chinDanishData.pos}')));