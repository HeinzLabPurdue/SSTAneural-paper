% function Fig4_compare_efficiency_difcor_c_t(saveFig, LatexDir)
function Fig4_compare_efficiency_difcor_c_t(saveFig, LatexDir)

if ~exist('saveFig', 'var')
    saveFig= 0;
end
if ~exist('LatexDir', 'var')
    LatexDir= ['figures' filesep];
end

dirStruct.latexDir= LatexDir;
dirStruct.loading_dir= ['data' filesep 'DanishData' filesep];


figSize_cm= [15 5 13.2 8];
figure_prop_name = {'PaperPositionMode','units','Position', 'Renderer'};
figure_prop_val =  { 'auto'            ,'centimeters', figSize_cm, 'painters'};  % [Xcorner Ycorner Xwidth Ywidth]
figure(1);
clf;
set(gcf,figure_prop_name,figure_prop_val);

[sig, fsOrg]= audioread(['stimuli' filesep 'FLN_Stim_S_P.wav']);
fsSig= 20e3;
sig= helper.gen_resample(sig, fsOrg, fsSig);
sig= helper.gen_rescale(sig, 65);



nh_hinIDs=[321 322 325 338 341 343 346 347 354 355 373 379];

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
allChinSpikeData= allChinSpikeData(strcmp({allChinSpikeData.noise}, 'SSN'));
allChinSpikeData= allChinSpikeData(ismember([allChinSpikeData.chinID], nh_hinIDs));
all_chin_track_unit_spl= [ [allChinSpikeData.chinID]', [allChinSpikeData.track]', [allChinSpikeData.unit]', [allChinSpikeData.SPL]'];
uniq_chin_track_unit_spl= unique(all_chin_track_unit_spl, 'rows');

%% check units with 75 reps
final_data= struct('spikes_pos', [], 'spikes_neg', [], 'tc', struct([]), 'CF_Hz', [], 'SR', []);
new_ind= 0;
for uniqVar= 1:size(uniq_chin_track_unit_spl, 1)
    cur_unit_inds= find(ismember(all_chin_track_unit_spl, uniq_chin_track_unit_spl(uniqVar, :), 'rows'));
    
    temp_spike= struct('spikes_pos', [], 'spikes_neg', []);
    for condVar= 1:length(cur_unit_inds)
        cur_cond_ind= cur_unit_inds(condVar);
        temp_spike.spikes_pos= [temp_spike.spikes_pos; allChinSpikeData(cur_cond_ind).SpikeTrains{1,1}];
        temp_spike.spikes_neg= [temp_spike.spikes_neg; allChinSpikeData(cur_cond_ind).SpikeTrains{1,2}];
    end
    if numel(temp_spike.spikes_pos)>=75 && allChinSpikeData(cur_cond_ind).CF_Hz<2e3 && allChinSpikeData(cur_cond_ind).CF_Hz>0.3e3
        new_ind= new_ind+1;
        final_data(new_ind).spikes_pos= temp_spike.spikes_pos;
        final_data(new_ind).spikes_neg= temp_spike.spikes_neg;
        final_data(new_ind).tc= allChinSpikeData(cur_cond_ind).TC;
        final_data(new_ind).CF_Hz= allChinSpikeData(cur_cond_ind).CF_Hz;
        final_data(new_ind).SR= allChinSpikeData(cur_cond_ind).SR;
    end
end

%% Now compare fractional power for difcor and c(t)
anl.FixedDelay= 5e-3; % This is approximate
anl.tStart= 225e-3; anl.tEnd= 325e-3;
anl.dur= anl.tEnd - anl.tStart;
anl.BinRes= 50e-6;
anl.tBinEdges= anl.tStart:anl.BinRes:anl.tEnd;
anl.t_uRate= (anl.tBinEdges(1:end-1)+anl.tBinEdges(2:end))/2;
anl.nw= 1.5;
anl.yrange= 50;
anl.nBoots= 12;
anl.doPlot= false; %true  | false;
anl.nRepsinBoot= 25;
anl.ValidFreqs= 585; %[580 590];
anl.nfft= 2^(nextpow2(2*length(anl.tBinEdges)));

rng('default');

unitPerform_difcor= nan(length(final_data), anl.nBoots);
unitPerform_comp_dpss= nan(length(final_data), anl.nBoots);
for unitVar= 1:length(final_data)
    curSpikeTrain_pos= final_data(unitVar).spikes_pos;
    curSpikeTrain_neg= final_data(unitVar).spikes_neg;
    
    unit_frac_dft= nan(1, anl.nBoots);
    unit_frac_dpss= nan(1, anl.nBoots);
    
    for bootVar= 1:anl.nBoots
        %     for bootVar= 1%:anl.nBoots
        bootInds_pos= randsample(numel(curSpikeTrain_pos), anl.nRepsinBoot);
        bootInds_neg= randsample(numel(curSpikeTrain_neg), anl.nRepsinBoot);
        
        boot_curSpikeTrain_pos= cellfun(@(x,y,z) x(x>y & x<z), curSpikeTrain_pos(bootInds_pos), repmat({anl.tStart+anl.FixedDelay}, numel(bootInds_pos), 1), ...
            repmat({anl.tEnd+anl.FixedDelay}, numel(bootInds_pos), 1), 'UniformOutput', false);
        boot_curSpikeTrain_neg= cellfun(@(x,y,z) x(x>y & x<z), curSpikeTrain_neg(bootInds_neg), repmat({anl.tStart+anl.FixedDelay}, numel(bootInds_pos), 1), ...
            repmat({anl.tEnd+anl.FixedDelay}, numel(bootInds_pos), 1), 'UniformOutput', false);
        
        uRate_pos= histcounts(cell2mat(boot_curSpikeTrain_pos), anl.tBinEdges);
        uRate_neg= histcounts(cell2mat(boot_curSpikeTrain_neg), anl.tBinEdges);
        uRate_comp= (uRate_pos-uRate_neg)/2;
        
        [NSACp,NSACdelays,AVGrate,TOTALspikes]= helper.SACfull_m(boot_curSpikeTrain_pos,anl.BinRes,anl.dur);
        NSACn= helper.SACfull_m(boot_curSpikeTrain_neg,anl.BinRes,anl.dur);
        NSCC= helper.SCCfull_m({boot_curSpikeTrain_pos, boot_curSpikeTrain_neg},anl.BinRes,anl.dur);
        NSAC= (NSACp+NSACn)/2;
        
        difcor= (NSCC-NSAC);
        
        [Pdft, Fdft]= helper.plot_dft(difcor, 1/anl.BinRes, 'yrange', anl.yrange, 'plot', anl.doPlot, 'nfft', anl.nfft);
        valid_dft_freqs= dsearchn(Fdft(:), anl.ValidFreqs);
        unit_frac_dft(bootVar)= db((sum(db2mag(Pdft(valid_dft_freqs))))/(sum(db2mag(Pdft)))); % db2mag because fft(acf(x))=fft(x)^2
        
        [Pxx, Fxx]= helper.plot_dpss_psd(uRate_comp, 1/anl.BinRes, 'nw', anl.nw, 'yrange', anl.yrange, 'plot', anl.doPlot, 'nfft', anl.nfft);
        
        valid_dpss_freqs= dsearchn(Fxx(:), anl.ValidFreqs);
        unit_frac_dpss(bootVar)= db((sum(db2pow(Pxx(valid_dpss_freqs))))/(sum(db2pow(Pxx))));
        
end
    unitPerform_difcor(unitVar, :)= unit_frac_dft;
    unitPerform_comp_dpss(unitVar, :)= unit_frac_dpss;
end

%%
F0_freq= 100;
F1_freq= 630;

plt.mrkSize= 6;
plt.mrkSize2= 8;
plt.lw3= 1.5;
plt.lw2_5= 1.25;
plt.lw2= 1;
plt.tick_len= [.025 .025];

figure(1);
clf;
xtick_vals_unit= 1:size(unitPerform_difcor,1);
xtick_labs_unit= cellfun(@(x) num2str(x), num2cell(xtick_vals_unit), 'UniformOutput', false);

xtick_vals_freq= [100 500 1e3];
xtick_labs_freq= cellfun(@(x) num2str(x), num2cell(xtick_vals_freq), 'UniformOutput', false);

ax(1)= subplot(221);
hold on;
[Pdft_sig, Fdft_sig, dftHan]= helper.plot_dft(sig(round(anl.tStart*fsSig):round(anl.tEnd*fsSig)), fsSig, 'yscale', 'dbspl', 'yrange', 40);
dft_close_ind= dsearchn(Fdft_sig(:), anl.ValidFreqs);
plot(Fdft_sig(dft_close_ind), 7.5+Pdft_sig(dft_close_ind), 'v', 'color', helper.get_color('g'), 'LineWidth', plt.lw3, 'markersize', plt.mrkSize);
text(F0_freq, max(Pdft_sig)+8, 'F_0', 'HorizontalAlignment', 'center');
text(F1_freq, max(Pdft_sig)+8, 'F_1', 'HorizontalAlignment', 'center');
xlabel('');
ylim([25 70])
set(gca, 'xscale', 'log', 'XTick', xtick_vals_freq, 'XTickLabel', xtick_labs_freq);
set(dftHan, 'linew', plt.lw2_5, 'color', 'k');
ylabHan(1)= ylabel('DFT (dB SPL)', 'Units', 'normalized');

ax(2)= subplot(223);
hold on;
lg(2)= plot(Fdft(2:end), db((db2mag(Pdft(2:end)))/(sum(db2mag(Pdft)))), 'LineWidth', plt.lw2); % 2:end because don't plot dc
lg(1)= plot(Fxx, db((db2pow(Pxx))/(sum(db2pow(Pxx)))), 'LineWidth', plt.lw2);
plot(Fdft(valid_dft_freqs), 10+db((db2mag(Pdft(valid_dft_freqs)))/(sum(db2mag(Pdft)))), 'v', 'color', helper.get_color('g'), 'LineWidth', plt.lw3, 'markersize', plt.mrkSize);
% plot(Fxx(valid_dpss_freqs), db((db2pow(Pxx(valid_dpss_freqs)))/(sum(db2pow(Pxx)))), '^', 'LineWidth', 3);
set(gca, 'xscale', 'log', 'XTick', xtick_vals_freq, 'XTickLabel', xtick_labs_freq);
set(lg(2), 'color', helper.get_color('r'));
set(lg(1), 'color', helper.get_color('b'));
% xlim([400  1e3])
xlim([80 1.1e3])
ylim([-80 -30])
%legend(lg, '$c(t)$', '$difcor$', 'box', 'off', 'interpreter', 'latex', 'Location', 'northwest');
xlabel('Frequency (Hz)')
ylabHan(2)= ylabel('PSD (dB/Hz)', 'Units', 'normalized');

ylabHan(1).Position(1)= mean([ylabHan(1).Position(1), ylabHan(2).Position(1)]);
ylabHan(2).Position(1)= mean([ylabHan(1).Position(1), ylabHan(2).Position(1)]);

ax(3)= subplot(222);
hold on;
lHan(1)= errorbar(1:length(final_data), nanmean(unitPerform_comp_dpss, 2), nanstd(unitPerform_comp_dpss, [], 2), 'LineWidth', plt.lw2, 'Color', helper.get_color('b'), 'LineStyle', '-');
lHan(2)= errorbar(1:length(final_data), nanmean(unitPerform_difcor, 2), nanstd(unitPerform_difcor, [], 2), 'LineWidth', plt.lw2, 'Color', helper.get_color('r'), 'LineStyle', '-');
set(gca, 'YColor', 'k');
ylabel('Power_{Frac}(6F_0) in dB');
[lgHan, icons]= legend(lg, 'D(f)', 'difcor', 'box', 'off', 'interpreter', 'tex', 'Location', 'northwest');
set(gca, 'XTick', xtick_vals_unit, 'XTickLabel', xtick_labs_unit);

icons(3).XData= mean(icons(3).XData) + [0 +.25];
icons(5).XData= mean(icons(5).XData) + [0 +.25];
lgHan.Position(1:2)= [.05 .36];

% lHan(3)= plot(nan, nan, 'dk', 'LineWidth', plt.lw3, 'MarkerSize', plt.mrkSize);

ax(4)= subplot(224);
hold on;
plot(1:length(final_data), (nanstd(unitPerform_difcor, [], 2) ./ nanstd(unitPerform_comp_dpss, [], 2)).^2, 'dk', 'LineWidth', plt.lw3, 'MarkerSize', plt.mrkSize2);
plot([0 1+length(final_data)], [1 1], 'color', 1.5*helper.get_color('gray'), 'LineWidth', 3, 'MarkerSize', plt.mrkSize2, 'LineStyle', '--');
set(gca, 'YColor', 'k');
ylabel('\sigma^2_{difcor}/\sigma^2_{D(f)}');
ylim([.5 2.5])

set(gca, 'XTick', xtick_vals_unit, 'XTickLabel', xtick_labs_unit);
% set(gcf, 'units', 'inches', 'position', [21 1 12 6.5]);
xlabel('Unit Number');
% legend(lHan, '$c(t)$', '$difcor$', '$var$-$ratio$', 'box', 'off', 'interpreter', 'latex', 'Location', 'southeast', 'FontSize', 14)

set(findall(gcf,'-property','FontSize'),'FontSize', 9);
set(findall(gcf,'-property','TickLength'),'TickLength', plt.tick_len, 'units', 'normalized');
helper.add_subplot_letter_txp(2, 2, 11, .05, 1.05);
linkaxes(ax(3:4), 'x')
xlim(ax(4), [min(xtick_vals_unit)-.5 max(xtick_vals_unit)+.5])
linkaxes(ax(1:2), 'x')
xlim(ax(1), [80 1.1e3])
%% define new axes for AB
Xcorner_AB= .08;
Xwidth_AB= .395;
Xshift_horz= .12;

Ycorner_AB= .11;
Ywidth_AB= .375;
Yshift_AB= .1;

% B
set(ax(2),'Position',[Xcorner_AB Ycorner_AB Xwidth_AB Ywidth_AB])
drawnow
% A
set(ax(1),'Position',[Xcorner_AB Ycorner_AB+Ywidth_AB+Yshift_AB Xwidth_AB Ywidth_AB])
drawnow
% define new axes for AB
Xwidth_CD=Xwidth_AB;
Ywidth_CD= Ywidth_AB;
Xcorner_CD= Xcorner_AB+Xwidth_AB+Xshift_horz;
Yshift_CD= Yshift_AB;
Ycorner_CD= Ycorner_AB;
% Yshift=0.06;

% D
set(ax(4),'Position',[Xcorner_CD Ycorner_CD Xwidth_CD Ywidth_CD])
drawnow
% C
set(ax(3),'Position',[Xcorner_CD Ycorner_CD+Ywidth_CD+Yshift_CD Xwidth_CD Ywidth_CD])
drawnow

% fName= sprintf('%sSingleFreqImprovement', dirStruct.latexDir);
if saveFig
    fName= sprintf('%sFig4', dirStruct.latexDir);
    saveas(gcf, fName, 'epsc');
end