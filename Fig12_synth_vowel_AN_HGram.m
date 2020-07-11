% function Fig12_synth_vowel_AN_HGram(saveFig, LatexDir)
function Fig12_synth_vowel_AN_HGram(saveFig, LatexDir)

if ~exist('saveFig', 'var')
    saveFig= 0;
end
if ~exist('LatexDir', 'var')
    LatexDir= ['figures' filesep];
end



% Init params
figHan.time= 1;
figHan.psd= 2;

anl.stimDuration= 188e-3; % stim ends at 188 ms
anl.stimPeriod= 350e-3; % stim rep rate (stim on + stim off)
anl.fixed_delay= 7.5e-3; % Approx 8 ms delay to AN fiber responses
anl.fs= 20e3; % 1/(bin resolution) for uRate/histcount
anl.tEdge_hist= 0:1/anl.fs:anl.stimPeriod;
anl.tPlot= (anl.tEdge_hist(1:end-1) + anl.tEdge_hist(2:end))/2;
anl.tStart= 38e-3;
anl.tEnd= anl.stimDuration; % anl.stimDuration || (38e-3+64e-3);
anl.tRange= anl.tEnd-anl.tStart;
anl.tResp= anl.tPlot>anl.fixed_delay & anl.tPlot<(anl.stimDuration+anl.fixed_delay);
anl.tMask= anl.tPlot>anl.tStart & anl.tPlot<anl.tEnd;
anl.kin_trajectory.f0= linspace(100, 120, round(anl.stimDuration*anl.fs));
anl.kin_trajectory.f1= linspace(630, 570, round(anl.stimDuration*anl.fs));
anl.kin_trajectory.f2= linspace(1200, 1500, round(anl.stimDuration*anl.fs));
anl.stat_params.f0= 100;
anl.stat_params.f1= 600;
anl.stat_params.f2= 1200;
anl.feature= 'RAW';
anl.nHarmonics= 15;

if length(anl.kin_trajectory.f0)~=sum(anl.tPlot<188e-3)
    error('Should be equal in length');
end

% Load saved data
DirStruct.INdata= ['data' filesep];
ChinID= 374;
SynCapData= load(sprintf('%sQ%d_allSyncCap.mat', DirStruct.INdata, ChinID));
SynCapData=SynCapData.SynCapData;

kin_Image= cell(length(SynCapData), 1);
%% Main loop to iterate through all units

d_lp = designfilt('lowpassiir','FilterOrder', 3, ...
    'HalfPowerFrequency', (1.5/anl.tRange)/(anl.fs/2), 'DesignMethod','butter');

parfor unitVar= 1:length(SynCapData)
    
    %% define what features to look at for stationary and kinematic vowels
    cur_unit_data= SynCapData(unitVar);
    if ~isempty(cur_unit_data.stat.RAW.pos)
        
        maxVal= -inf;
        
        %% kinematic
        temp_kin_uRate_pos= histcounts(cell2mat(cur_unit_data.kin.(anl.feature).pos'), anl.tEdge_hist);
        temp_kin_uRate_neg= histcounts(cell2mat(cur_unit_data.kin.(anl.feature).neg'), anl.tEdge_hist);
        temp_kin_uRate_compound= (temp_kin_uRate_pos-temp_kin_uRate_neg)/2;
        temp_kin_uRate_HilbMag= abs(hilbert(temp_kin_uRate_compound));
        temp_kin_uRate_HilbPhi= temp_kin_uRate_compound ./ temp_kin_uRate_HilbMag;
        
        temp_kin_Image= nan(anl.nHarmonics, sum(anl.tResp));
        for harmVar= 1:anl.nHarmonics
            temp_kin_Image(harmVar, :)= helper.get_trajectory_signal(temp_kin_uRate_HilbPhi(anl.tResp), anl.fs, harmVar*anl.kin_trajectory.f0, d_lp);
        end
        kin_Image{unitVar}= temp_kin_Image;
    end
end

%%
all_cfs= [SynCapData.CF_Hz];

CF_boundary= [0 1e3 2.5e3 inf];
avg_kin_image= repmat({zeros(anl.nHarmonics, sum(anl.tMask))}, 3, 1);

for unitVar= 1:length(SynCapData)
    cf_Hz_cur= all_cfs(unitVar);
    cf_ind= cf_Hz_cur>CF_boundary(1:end-1) & cf_Hz_cur<CF_boundary(2:end);
    if ~isempty(kin_Image{unitVar})
        avg_kin_image{cf_ind}= avg_kin_image{cf_ind}+ kin_Image{unitVar}(:, anl.tMask);
    end
end


%%
figSize_cm= [5 5 13.2 6]; % [Xcorner Ycorner Xwidth Ywidth]
figHan.HGram= 11;
figure_prop_name = {'PaperPositionMode','units','Position'};
figure_prop_val =  { 'auto'            ,'centimeters', figSize_cm};  % [Xcorner Ycorner Xwidth Ywidth]
figure(figHan.HGram);
clf;
set(figHan.HGram, figure_prop_name, figure_prop_val);


lw= 1.5;
y_harmonics= (1:anl.nHarmonics);
x_tSamples= (anl.tPlot(anl.tMask))*1e3;
clf;
title_str= {'A. Pooled low CFs', 'B. Pooled medium CFs'};
txtHan= nan(2, 1);
sp_ax= nan(2, 1);
for cf_range_var= 1:2
    sp_ax(cf_range_var)= subplot(1, 2, cf_range_var);
    imagesc(x_tSamples, y_harmonics, avg_kin_image{cf_range_var})
    hold on;
    lineHan= line(x_tSamples, anl.kin_trajectory.f1(anl.tMask)./anl.kin_trajectory.f0(anl.tMask), 'LineWidth', lw);
    set(lineHan, 'color', helper.get_color('prp'), 'LineStyle', '-');
    lineHan= line(x_tSamples, ones(size(x_tSamples)), 'LineWidth', lw);
    set(lineHan, 'color', helper.get_color('k'), 'LineStyle', '-');
    if cf_range_var==2
        lineHan= line(x_tSamples, anl.kin_trajectory.f2(anl.tMask)./anl.kin_trajectory.f0(anl.tMask), 'LineWidth', lw);
        set(lineHan, 'color', helper.get_color('r'), 'LineStyle', '-');
    end
    
    set(gca, 'YDir', 'normal');
    txtHan(cf_range_var)= text(.05, 1.05, title_str{cf_range_var}, 'units', 'normalized');
end
subplot(1, 2, 1);
xlab_han= xlabel('Time (ms)');
xlab_han.Position([1 2])= [190 -.5];

subplot(1, 2, 1);
ylabel('Harmonic number');

linkaxes(sp_ax);
set(findall(gcf,'-property','FontSize'),'FontSize', 9);
set(txtHan,'FontSize', 11);

% Set axes placement/size
Xwidth=.42;
Xcorner=.075;
Xshift=.075;
Ywidth=.8;
Ycorner=.12;

% A
set(sp_ax(1),'Position',[Xcorner Ycorner Xwidth Ywidth])
set(txtHan(1),'pos',get(txtHan(1),'pos')+[0 0.01 0])
drawnow
% B
set(sp_ax(2),'Position',[Xcorner+Xwidth+Xshift Ycorner Xwidth Ywidth])
set(txtHan(2),'pos',get(txtHan(2),'pos')+[0 0.01 0])
drawnow

if saveFig
    saveas(gcf, [LatexDir 'Fig12'], 'epsc');
end