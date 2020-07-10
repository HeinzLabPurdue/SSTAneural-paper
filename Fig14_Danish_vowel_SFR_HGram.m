% function Fig14_Danish_vowel_SFR_HGram(saveFig, LatexDir)
function Fig14_Danish_vowel_SFR_HGram(saveFig, LatexDir)

if ~exist('saveFig', 'var')
    saveFig= 0;
end
if ~exist('LatexDir', 'var')
    LatexDir= ['figures' filesep];
end

file2load= ['data' filesep 'SP-2019_06_09-Q374_SFR_NH' filesep 'a0008_FFR_SNRenvSSN_Stim_S_P_atn10.mat'];

hilbPhi1_comp0= 0;

%%
data= load(file2load);
data= data.data;

anl.fs= 10e3;
fs_data= data.Stimuli.RPsamprate_Hz;
gain= 20e3;

x_pos= data.AD_Data.AD_Avg_PO_V{1} * (1e6/gain); % uV
x_neg= data.AD_Data.AD_Avg_NP_V{1} * (1e6/gain); % uV


plt.tWindow= 64e-3;
plt.fracOVlap= .98;
plt.nfft= 2^(1+nextpow2(anl.fs*plt.tWindow));
plt.doPlot= 0;
plt.lw3= 1.5;
plt.lw2= 1.0;

x_tfs= (x_pos-x_neg)/2;
x_tfs= helper.gen_resample(x_tfs, fs_data, anl.fs);
x_hil= hilbert(x_tfs);
x_hil_mag= abs(x_hil);
if hilbPhi1_comp0
    x_TFS= x_hil ./ x_hil_mag;
else
    x_TFS= x_tfs;
end
t_ffr= (1:length(x_tfs))/anl.fs;


x_env= (x_pos+x_neg)/2;
x_env= helper.gen_resample(x_env, fs_data, anl.fs);

%%
[sig, fs_sig]= audioread(['stimuli' filesep 'FLN_Stim_S_P.wav']);
sig= helper.gen_resample(sig, fs_sig, anl.fs);

anl.stimDuration= numel(sig)/anl.fs;
anl.fixed_delay= 7.5e-3; % Approx 8 ms delay to AN fiber responses
anl.nHarmonics= 20;
anl.tPlot= (1:length(x_TFS))/anl.fs;
anl.tStart= 250e-3; %200e-3;  | 860e-3
anl.tEnd= .66; %750e-3; %anl.stimDuration; % anl.stimDuration || 1250e-3
anl.tRange= anl.tEnd-anl.tStart;
anl.tResp= anl.tPlot>anl.fixed_delay & anl.tPlot<(anl.stimDuration+anl.fixed_delay); % This takes care of fixed delay due to inv FIR filtering
anl.tMask= anl.tPlot>anl.tStart & anl.tPlot<anl.tEnd;

d_lp = designfilt('lowpassiir','FilterOrder', 2, ...
    'HalfPowerFrequency', (10/anl.tRange)/(anl.fs/2), 'DesignMethod','butter');

%%
danish.voiced_boundaries= helper.find_voicing_boundaries(sig, anl.fs, 0, .13);

temp_f0= load(['data' filesep 'danish_pitch.mat']);
temp_f0= temp_f0.pitch_data;
danish.voiced_inds= any(anl.tPlot(anl.tResp)>danish.voiced_boundaries(:,1) & anl.tPlot(anl.tResp)<danish.voiced_boundaries(:,2), 1);
danish.trajectory.f0= zeros(size(anl.tResp(anl.tResp)));
danish.trajectory.f0(danish.voiced_inds)= interp1([temp_f0.time], [temp_f0.est], anl.tPlot(danish.voiced_inds), 'pchip');
danish.voiced_inds= danish.voiced_inds & danish.trajectory.f0>95 & danish.trajectory.f0<150; % make talker - looking at estiamtes, This removes the edge estimates.
danish.trajectory.f0(~danish.voiced_inds)= 0;

temp_formants= load(['data' filesep 'danish_formant.mat']);
temp_formants= temp_formants.formant_data;
danish.trajectory.f1= interp1([temp_formants.time], [temp_formants.f1], anl.tPlot, 'pchip');
danish.trajectory.f1(~danish.voiced_inds)= nan;
danish.trajectory.f1(round(anl.fs*anl.stimDuration):end)= nan;
danish.trajectory.f2= interp1([temp_formants.time], [temp_formants.f2], anl.tPlot, 'pchip');
danish.trajectory.f2(~danish.voiced_inds)= nan;
danish.trajectory.f2(round(anl.fs*anl.stimDuration):end)= nan;
danish.trajectory.f3= interp1([temp_formants.time], [temp_formants.f3], anl.tPlot, 'pchip');
danish.trajectory.f3(~danish.voiced_inds)= nan;
danish.trajectory.f3(round(anl.fs*anl.stimDuration):end)= nan;


temp_danish_Image= nan(anl.nHarmonics, sum(anl.tResp));
for harmVar= 1:anl.nHarmonics
    temp_danish_Image(harmVar, :)= helper.get_trajectory_signal(x_TFS(anl.tResp), anl.fs, harmVar*danish.trajectory.f0, d_lp);
end

x_tSamples= (anl.tPlot(anl.tMask));
plt.spect_xshift= min(x_tSamples);
y_harmonics= (1:anl.nHarmonics);

%%
figSize_cm= [5 5 13.2 10]; % [Xcorner Ycorner Xwidth Ywidth]
figHan= 1;
figure_prop_name = {'PaperPositionMode','units','Position'};
figure_prop_val =  { 'auto'            ,'centimeters', figSize_cm};  % [Xcorner Ycorner Xwidth Ywidth]
figure(figHan);
clf;
set(figHan, figure_prop_name, figure_prop_val);

xtick_vals= 300:100:600;
SPhan(1)= subplot(211);
helper.plot_spectrogram(x_TFS, anl.fs, plt.tWindow, plt.fracOVlap, plt.nfft, plt.doPlot);
plt.ColBar_Spec= max(get(colorbar, 'Limits'));
if ~hilbPhi1_comp0
    caxis([plt.ColBar_Spec-15 plt.ColBar_Spec]);
else
    caxis([plt.ColBar_Spec-14 plt.ColBar_Spec]);
end
colorbar off;
txtHan(1)= text(.05, 1.1, 'A. FFR Spectrogram', 'units', 'normalized');
ylim([0 2.2]);
ylHan(1)= ylabel('Frequency (kHz)', 'Units', 'normalized');

hold on;
plot(x_tSamples*1e3, danish.trajectory.f1(anl.tMask)/1e3, 'color', helper.get_color('prp'),'LineWidth', plt.lw3);
plot(x_tSamples*1e3, danish.trajectory.f2(anl.tMask)/1e3, 'color', helper.get_color('r'), 'LineWidth', plt.lw3);

yyaxis right;
plot(x_tSamples*1e3, danish.trajectory.f0(anl.tMask), 'k-', 'LineWidth', plt.lw3);
set(gca, 'YColor', 'k', 'XTick', xtick_vals);
ylabel('F_0 Frequency (Hz)');
xlabel('');

SPhan(2)= subplot(212);
imagesc(x_tSamples*1e3, y_harmonics, db(temp_danish_Image(:, anl.tMask)))
set(gca, 'YDir', 'normal', 'XTick', xtick_vals);
hold on;
plot((x_tSamples*1e3), danish.trajectory.f1(anl.tMask)./danish.trajectory.f0(anl.tMask), 'color', helper.get_color('prp'), 'LineWidth', plt.lw3);
plot((x_tSamples*1e3), danish.trajectory.f2(anl.tMask)./danish.trajectory.f0(anl.tMask), 'color', helper.get_color('r'), 'LineWidth', plt.lw3);
plot((x_tSamples*1e3), danish.trajectory.f0(anl.tMask)./danish.trajectory.f0(anl.tMask), 'k-', 'LineWidth', plt.lw3);
ylHan(2)= ylabel('Harmonic Number', 'Units', 'normalized');
txtHan(2)= text(.05, 1.1, 'B. FFR Harmonicgram', 'units', 'normalized');
xlabel('Time (ms)');

plt.ColBar_Harm= max(get(colorbar, 'Limits'));
caxis([plt.ColBar_Harm-21 plt.ColBar_Harm]);
colorbar off;

linkaxes(SPhan, 'x');
xlim([anl.tStart anl.tEnd]*1e3)
ylHan(1).Position(1)= min([ylHan(1).Position(1) ylHan(2).Position(1)]);
ylHan(2).Position(1)= min([ylHan(1).Position(1) ylHan(2).Position(1)]);

% Set axes placement/size
Xcorner=.085;
Xwidth=.81;
Ywidth=.36;
Yshift=.115;
Ycorner=.1;

% A
set(SPhan(2),'Position',[Xcorner Ycorner Xwidth Ywidth])
drawnow
% B
set(SPhan(1),'Position',[Xcorner Ycorner+Ywidth+Yshift Xwidth Ywidth])
drawnow

set(findall(gcf,'-property','FontSize'),'FontSize', 9);
set(txtHan,'FontSize', 11);

if hilbPhi1_comp0
    fName= 'SigS3';
else
    fName= 'Fig14';
end
if saveFig
    print([LatexDir fName], '-dpng',  '-r600');
end