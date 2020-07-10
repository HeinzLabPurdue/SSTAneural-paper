% function Fig13_Danish_vowel_AN_HGram(saveFig, LatexDir)
function Fig13_Danish_vowel_AN_HGram(saveFig, LatexDir)

if ~exist('saveFig', 'var')
    saveFig= 0;
end
if ~exist('LatexDir', 'var')
    LatexDir= ['figures' filesep];
end

% Load saved data
DirStruct.INdata= ['data' filesep];
DirStruct.latex= LatexDir;

anl.stimDuration= 1300e-3; % stim ends at 188 ms
anl.stimPeriod= 1800e-3; % stim rep rate (stim on + stim off)
anl.fixed_delay= 7.5e-3; % Approx 8 ms delay to AN fiber responses
anl.fs= 20e3; % 1/(bin resolution) for uRate/histcount
anl.tEdge_hist= 0:1/anl.fs:anl.stimPeriod;
anl.tBinCenter= (anl.tEdge_hist(1:end-1)+anl.tEdge_hist(2:end))/2;
anl.tPlot= (anl.tEdge_hist(1:end-1) + anl.tEdge_hist(2:end))/2;
anl.tStart= 200e-3; %200e-3;  | 860e-3
anl.tEnd= 750e-3; %750e-3; %anl.stimDuration; % anl.stimDuration || 1250e-3
anl.tRange= anl.tEnd-anl.tStart;
anl.tResp= anl.tPlot>anl.fixed_delay & anl.tPlot<(anl.stimDuration+anl.fixed_delay); % This takes care of fixed delay due to inv FIR filtering
anl.tMask= anl.tPlot>anl.tStart & anl.tPlot<anl.tEnd;

[sig, fsSig]= audioread(['stimuli' filesep 'FLN_Stim_S_P.wav']);
danish.voiced_boundaries= helper.find_voicing_boundaries(sig, fsSig, 0, .13);

temp_f0= load(['data' filesep 'danish_pitch.mat']);
temp_f0= temp_f0.pitch_data;
danish.voiced_inds= any(anl.tBinCenter(anl.tResp)>danish.voiced_boundaries(:,1) & anl.tBinCenter(anl.tResp)<danish.voiced_boundaries(:,2), 1);
danish.trajectory.f0= zeros(size(anl.tBinCenter(anl.tResp)));
danish.trajectory.f0(danish.voiced_inds)= interp1([temp_f0.time], [temp_f0.est], anl.tBinCenter(danish.voiced_inds), 'pchip');
danish.voiced_inds= danish.voiced_inds & danish.trajectory.f0>95 & danish.trajectory.f0<150; % make talker - looking at estiamtes, This removes the edge estimates.
danish.trajectory.f0(~danish.voiced_inds)= 0;


temp_formants= load(['data' filesep 'danish_formant.mat']);
temp_formants= temp_formants.formant_data;
danish.trajectory.f1= interp1([temp_formants.time], [temp_formants.f1], anl.tBinCenter(anl.tResp), 'pchip');
danish.trajectory.f1(~danish.voiced_inds)= nan;
danish.trajectory.f2= interp1([temp_formants.time], [temp_formants.f2], anl.tBinCenter(anl.tResp), 'pchip');
danish.trajectory.f2(~danish.voiced_inds)= nan;
danish.trajectory.f3= interp1([temp_formants.time], [temp_formants.f3], anl.tBinCenter(anl.tResp), 'pchip');
danish.trajectory.f3(~danish.voiced_inds)= nan;

anl.feature= 'RAW';
anl.nHarmonics= 35;
anl.nHarmonics2Plot= 20;

if length(danish.trajectory.f0)~=sum(anl.tPlot(anl.tResp)<1800e-3)
    error('Should be equal in length');
end

ChinID= 374;
SynCapData= load(sprintf('%sQ%d_allSyncCap.mat', DirStruct.INdata, ChinID));
SynCapData=SynCapData.SynCapData;


allUnit_vowelPSDcell= cell(length(SynCapData), 1);

danish_Image= cell(length(SynCapData), 1);
danish_F1= nan(length(SynCapData), sum(anl.tResp));
danish_F2= nan(length(SynCapData), sum(anl.tResp));
danish_NF= nan(length(SynCapData), sum(anl.tResp));

%% Main loop to iterate through all units

d_lp_f0= designfilt('lowpassiir','FilterOrder', 3, ...
    'HalfPowerFrequency', (2.5/anl.tRange)/(anl.fs/2), 'DesignMethod','butter');
d_lp_formant= designfilt('lowpassiir','FilterOrder', 3, ...
    'HalfPowerFrequency', 200/(anl.fs/2), 'DesignMethod','butter');

uRate_diffs= nan(length(SynCapData), length(anl.tBinCenter));

parfor unitVar= 1:length(SynCapData)
    
    %% define what features to look at for danish
    cur_unit_data= SynCapData(unitVar);
    if ~isempty(cur_unit_data.kin.RAW.pos)
        maxVal= -inf;
        
        %% kinematic
        temp_danish_uRate_pos= histcounts(cell2mat(cur_unit_data.danish.pos), anl.tEdge_hist);
        temp_danish_uRate_neg= histcounts(cell2mat(cur_unit_data.danish.neg), anl.tEdge_hist);
        temp_danish_uRate_compound= (temp_danish_uRate_pos-temp_danish_uRate_neg)/2;
        temp_danish_uRate_HilbMag= abs(hilbert(temp_danish_uRate_compound));
        temp_danish_uRate_HilbPhi= rms(temp_danish_uRate_compound) * temp_danish_uRate_compound ./ temp_danish_uRate_HilbMag;
        uRate_diffs(unitVar, :)= temp_danish_uRate_HilbPhi;
        
        temp_danish_Image= nan(anl.nHarmonics, sum(anl.tResp));
        for harmVar= 1:anl.nHarmonics
            temp_danish_Image(harmVar, :)= helper.get_trajectory_signal(temp_danish_uRate_HilbPhi(anl.tResp), anl.fs, harmVar*danish.trajectory.f0, d_lp_f0);
        end
        danish_Image{unitVar}= temp_danish_Image;
        
        [~, danish_F1(unitVar, :)]= helper.get_trajectory_signal(temp_danish_uRate_HilbPhi(anl.tResp), anl.fs, danish.trajectory.f1, d_lp_formant);
        [~, danish_F2(unitVar, :)]= helper.get_trajectory_signal(temp_danish_uRate_HilbPhi(anl.tResp), anl.fs, danish.trajectory.f2, d_lp_formant);
        [~, danish_NF(unitVar, :)]= helper.get_trajectory_signal(temp_danish_uRate_HilbPhi(anl.tResp), anl.fs, 3e3*ones(size(danish.trajectory.f2)), d_lp_formant);
    end
end

%%
all_cfs= [SynCapData.CF_Hz];

CF_boundary= [0 1e3 2.5e3 inf];
CF_labels= {'LF', 'MF', 'HF'};

avg_danish_image= repmat({zeros(anl.nHarmonics, sum(anl.tMask))}, 3, 1);
for unitVar= 1:length(SynCapData)
    cf_Hz_cur= all_cfs(unitVar);
    cf_ind= cf_Hz_cur>CF_boundary(1:end-1) & cf_Hz_cur<CF_boundary(2:end);
    if ~isempty(danish_Image{unitVar})
        avg_danish_image{cf_ind}= avg_danish_image{cf_ind}+ danish_Image{unitVar}(:, anl.tMask);
    end
end

%% compare harmonic plot versus spectrogram
useF1_0_or_F1harm_1= 1; % If 0 use 400 Hz BW around Fx, if 1 add three harmonics around Fx/F0.
plt.lw1= .85;
plt.lw2= 1.25;
plt.tick_len= [.025 .025];
y_harmonics= (1:anl.nHarmonics2Plot);
x_tSamples= (anl.tPlot(anl.tMask))*1e3;


figSize_cm= [15 5 13.2 10];
figure_prop_name = {'PaperPositionMode', 'units', 'Position', 'Renderer'};
figure_prop_val =  { 'auto', 'centimeters', figSize_cm, 'painters'};  % [Xcorner Ycorner Xwidth Ywidth]
figure(1);
clf;
set(gcf,figure_prop_name,figure_prop_val);


plt.ColBar= -inf;
plt.spect_default= 0;
if plt.spect_default
    plt.spect_xshift= min(x_tSamples);
else
    plt.spect_xshift= 0;
end

nSProws= 3;
nSPcols= 2;
sp_ax= nan(nSProws*nSPcols, 1);
titleStrs= {'A', 'B', 'C', 'D', 'E', 'F'};

danish.trajectory.f0(danish.trajectory.f0<50)= nan;
anl.f0_tMask= anl.tMask;
f0_traj_plot= danish.trajectory.f0;
f0_traj_plot(anl.tBinCenter<.25)= nan;
f0_traj_plot(anl.tBinCenter>.67)= nan;

harmMat= (1:anl.nHarmonics)' + zeros(size(danish.trajectory.f0(anl.tMask)));
f1_harmInds= danish.trajectory.f1(anl.tMask)./danish.trajectory.f0(anl.tMask);
f1_harmInds= [floor(f1_harmInds-.5); round(f1_harmInds); ceil(f1_harmInds+.5)]';
f2_harmInds= danish.trajectory.f2(anl.tMask)./danish.trajectory.f0(anl.tMask);
f2_harmInds= [floor(f2_harmInds-.5); round(f2_harmInds); ceil(f2_harmInds+.5)]';
nf_harmInds= [-1, 0, 1] + 30*ones(size(danish.trajectory.f0(anl.tMask)'));

nPoints= size(f1_harmInds, 1);
nHarm_near_cf= size(f1_harmInds, 2);

f1_harm_mask= zeros(size(avg_danish_image{1}));
f2_harm_mask= zeros(size(avg_danish_image{1}));
nf_harm_mask= zeros(size(avg_danish_image{1}));

temp_f1_subX= f1_harmInds(:);
temp_f1_subX(temp_f1_subX==0)= nan;
temp_f1_subY= repmat((1:nPoints)', nHarm_near_cf, 1);
temp_valid_f1inds= find(~(isnan(temp_f1_subX) | isnan(temp_f1_subY)));
temp_f1_inds= sub2ind(size(f1_harm_mask), temp_f1_subX(temp_valid_f1inds), temp_f1_subY(temp_valid_f1inds));
f1_harm_mask(temp_f1_inds)= 1;

temp_f2_subX= f2_harmInds(:);
temp_f2_subX(temp_f2_subX==0)= nan;
temp_f2_subY= repmat((1:nPoints)', nHarm_near_cf, 1);
temp_valid_f2inds= find(~(isnan(temp_f2_subX) | isnan(temp_f2_subY)));
temp_f2_inds= sub2ind(size(f2_harm_mask), temp_f2_subX(temp_valid_f2inds), temp_f2_subY(temp_valid_f2inds));
f2_harm_mask(temp_f2_inds)= 1;

temp_nf_subX= nf_harmInds(:);
temp_nf_subX(temp_nf_subX==0)= nan;
temp_nf_subY= repmat((1:nPoints)', nHarm_near_cf, 1);
temp_valid_nfinds= find(~(isnan(temp_nf_subX) | isnan(temp_nf_subY)));
temp_nf_inds= sub2ind(size(nf_harm_mask), temp_nf_subX(temp_valid_nfinds), temp_nf_subY(temp_valid_nfinds));
nf_harm_mask(temp_nf_inds)= 1;

for cf_range_var= 1:nSPcols% length(avg_danish_image)
    cur_cf_lims= CF_boundary(cf_range_var:cf_range_var+1);
    cur_cf_inds= all_cfs>cur_cf_lims(1) & all_cfs<cur_cf_lims(2);
    cur_avg_diffRate= nanmean(uRate_diffs(cur_cf_inds, :), 1);
    
    sp_ax(cf_range_var)= subplot(nSProws, nSPcols, cf_range_var);
    
    yyaxis right;
    lHan= line((x_tSamples), f0_traj_plot(anl.tMask), 'LineWidth', plt.lw2);
    set(lHan, 'color', helper.get_color('k'), 'LineStyle', '-');
    set(gca, 'YColor', 'k');
    if cf_range_var==1
        set(gca, 'YTickLabel', '');
    end
    set(gca, 'XTickLabel', '');
    
    yyaxis left;
    hold on;
    helper.plot_spectrogram(cur_avg_diffRate(anl.tMask), anl.fs, [], [], [], plt.spect_default, anl.tStart);
    plt.ColBar= max(plt.ColBar, max(get(colorbar, 'Limits')));
    
    lHan= line((x_tSamples), danish.trajectory.f1(anl.tMask)/1e3, 'LineWidth', plt.lw2);
    set(lHan, 'color', helper.get_color('prp'), 'LineStyle', '-');
    
    if cf_range_var==2
        lHan= line((x_tSamples), danish.trajectory.f2(anl.tMask)/1e3, 'LineWidth', plt.lw2);
        set(lHan, 'color', helper.get_color('r'), 'LineStyle', '-');
    end
    
    ylim([0 2.5]);
    text(.05, 1.1, titleStrs{cf_range_var}, 'units', 'normalized');
    set(gca, 'YColor', 'k');
    
    
    sp_ax(cf_range_var+nSPcols)= subplot(nSProws, nSPcols, cf_range_var+nSPcols);
    cur_cf_avgImage= avg_danish_image{cf_range_var};
    imagesc(x_tSamples, y_harmonics, cur_cf_avgImage(y_harmonics, :))
    hold on;
    lHan= line((x_tSamples - plt.spect_xshift + 32), ones(size(find(anl.tMask))), 'LineWidth', plt.lw2);
    set(lHan, 'color', helper.get_color('k'), 'LineStyle', '-');
    
    lHan= line(x_tSamples, danish.trajectory.f1(anl.tMask)./danish.trajectory.f0(anl.tMask), 'LineWidth', plt.lw2);
    set(lHan, 'color', helper.get_color('prp'), 'LineStyle', '-');
    
    if cf_range_var~=1
        lHan= line(x_tSamples, danish.trajectory.f2(anl.tMask)./danish.trajectory.f0(anl.tMask), 'LineWidth', plt.lw2);
        set(lHan, 'color', helper.get_color('r'), 'LineStyle', '-');
    end
    set(gca, 'YDir', 'normal');
    set(gca, 'XTickLabel', '');
    text(.05, 1.1, titleStrs{cf_range_var+nSPcols}, 'units', 'normalized');
    xlim([231 699]);
    
    sp_ax(cf_range_var+2*nSPcols)= subplot(nSProws, nSPcols, cf_range_var+2*nSPcols);
    if useF1_0_or_F1harm_1
        cur_F1_power= helper.trifilt(nanmean(cur_cf_avgImage .* (f1_harm_mask==1), 1), 5);
        cur_F2_power= helper.trifilt(nanmean(cur_cf_avgImage .* (f2_harm_mask==1), 1), 5);
        cur_NF_power= helper.trifilt(nanmean(cur_cf_avgImage .* (nf_harm_mask==1), 1), 5);
        hold on;
        formHan(2)= plot((x_tSamples - plt.spect_xshift), cur_F2_power, 'LineWidth', plt.lw2, 'color', helper.get_color('r'));
        formHan(1)= plot((x_tSamples - plt.spect_xshift), cur_F1_power, 'LineWidth', plt.lw1, 'color', helper.get_color('prp'));
        formHan(3)= plot((x_tSamples - plt.spect_xshift), cur_NF_power, 'LineWidth', plt.lw1, 'color', helper.get_color('gray'));
    else
        cur_F1_power= nanmean(danish_F1(cur_cf_inds, :), 1);
        cur_F2_power= nanmean(danish_F2(cur_cf_inds, :), 1);
        cur_NF_power= nanmean(danish_NF(cur_cf_inds, :), 1);
        hold on;
        formHan(2)= plot((x_tSamples - plt.spect_xshift), cur_F2_power(anl.tMask), 'LineWidth', plt.lw2, 'color', helper.get_color('r'));
        formHan(1)= plot((x_tSamples - plt.spect_xshift), cur_F1_power(anl.tMask), 'LineWidth', plt.lw1, 'color', helper.get_color('b'));
        formHan(3)= plot((x_tSamples - plt.spect_xshift), cur_NF_power(anl.tMask), 'LineWidth', plt.lw1, 'color', helper.get_color('gray'));
    end
    text(.05, 1.1, titleStrs{cf_range_var+2*nSPcols}, 'units', 'normalized');
end
[lgHan, icons]= legend(formHan, 'F_1', 'F_2', 'NF');
lgHan.Box= 'off';

lgHan.Position(1:2)= [0.78, 0.2];

icons(4).XData= mean(icons(4).XData) + [0.1 .3];
icons(4).LineWidth= 1;

icons(6).XData= mean(icons(6).XData) + [0.1 .3];
icons(6).YData= icons(6).YData + .12;
icons(2).Position(2)= mean(icons(6).YData);
icons(6).LineWidth= 1;

icons(8).XData= mean(icons(8).XData) + [0.1 .3];
icons(8).YData= icons(8).YData + .18;
icons(3).Position(2)= mean(icons(8).YData);
icons(8).LineWidth= 1;

for cf_range_var=1:nSPcols
    subplot(nSProws, nSPcols, cf_range_var) ;
    caxis([plt.ColBar-15 plt.ColBar]);
    colorbar off;
    if cf_range_var~=3
        ylabel('')
        xlabel('')
    end
end

subplot(nSProws, nSPcols, 1);
ylabHan(1)= ylabel('Frequency (kHz)', 'Units', 'normalized');

subplot(nSProws, nSPcols, 3);
ylabHan(2)= ylabel('Harmonic number', 'Units', 'normalized');


subplot(nSProws, nSPcols, 2);
yyaxis right;
ylabel('F_0 Frequency (Hz)');

subplot(nSProws, nSPcols, 5);
ylabHan(3)= ylabel('Formant Power', 'Units', 'normalized');
xlab_han= xlabel('Time (ms)');
xlab_han.Position([1 2])= [710 -1];

ylab_X= min([ylabHan(1).Position(1), ylabHan(2).Position(1), ylabHan(3).Position(1)]);
ylabHan(1).Position(1)= ylab_X;
ylabHan(2).Position(1)= ylab_X;
ylabHan(3).Position(1)= ylab_X;

set(findall(gcf,'-property','FontSize'),'FontSize', 9);
linkaxes(sp_ax, 'x');
xlim(sp_ax(1), [231 669]);
lgHan.FontSize= 12;

linkaxes(sp_ax([5 6]), 'y');
set(findall(gcf,'-property','TickLength'),'TickLength', plt.tick_len, 'units', 'normalized');

ttlHan(1)= title(sp_ax(1), 'Pooled low CFs', 'Units', 'normalized');
ttlHan(2)= title(sp_ax(2), 'Pooled high CFs', 'Units', 'normalized');
ttlHan(1).Position(1:2)= [.5 1.08];
ttlHan(2).Position(1:2)= [.5 1.08];

%% define new axes for AB
Xcorner_X= .075;
Xwidth_X= .375;
Xshift_X= .07;
Ycorner_X= .08;
Ywidth_X= .24;
Yshift_X= .07;

% E
set(sp_ax(5),'Position',[Xcorner_X, Ycorner_X, Xwidth_X Ywidth_X])
drawnow
% C
set(sp_ax(3),'Position',[Xcorner_X, Ycorner_X+Ywidth_X+Yshift_X, Xwidth_X, Ywidth_X])
drawnow
% A
set(sp_ax(1),'Position',[Xcorner_X, Ycorner_X+2*Ywidth_X+2*Yshift_X, Xwidth_X, Ywidth_X])
drawnow
% F
set(sp_ax(6),'Position',[Xcorner_X+Xwidth_X+Xshift_X, Ycorner_X, Xwidth_X, Ywidth_X])
drawnow
% D
set(sp_ax(4),'Position',[Xcorner_X+Xwidth_X+Xshift_X, Ycorner_X+Ywidth_X+Yshift_X, Xwidth_X, Ywidth_X])
drawnow
% B
set(sp_ax(2),'Position',[Xcorner_X+Xwidth_X+Xshift_X, Ycorner_X+2*Ywidth_X+2*Yshift_X, Xwidth_X, Ywidth_X])
drawnow

if saveFig
    fName_latex= [DirStruct.latex 'Fig13'];
%     fName_latex_tiff= [DirStruct.latex 'tiff/Fig13'];
    print(fName_latex, '-dpng',  '-r600');
%     print(fName_latex_tiff, '-dtiff',  '-r600');
end