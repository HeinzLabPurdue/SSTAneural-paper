% function Fig6_7_RD_FM_example(saveFig, LatexDir)
function Fig6_7_RD_FM_example(saveFig, LatexDir)

if ~exist('saveFig', 'var')
    saveFig= 0;
end
if ~exist('LatexDir', 'var')
    LatexDir= ['figures' filesep];
end

modFreq= 20; %[10 25 50 100];

FigHan.ENV= 1;
FigHan.TFS= 2;

plt.ylim.vIHC= [-150 -80];
plt.ylim.SynOut= [-45 25];
plt.xlim= [50 7.01e3];
plt.tick_len1= [.04 .04];
plt.tick_len2= [.012 .012];
plt.xtick= [100 1e3 5e3];
plt.fSize= 9;
plt.titleFsize= 20;

StimParams.fs= 100e3;
StimParams.fm= modFreq;
StimParams.modDepth= 1;
StimParams.dur= 1;
StimParams.phi_m= [];
StimParams.dBreThresh= 40;
StimParams.rlf_dur= .4;
StimParams.dBSPL= 30;
StimParams.DrivenRateTarget= 130;

% model parameters
ANparams.spont = 70;   % spontaneous firing rate
ANparams.tabs= 0.6e-3; % Absolute refractory period
ANparams.trel= 0.6e-3; % Baseline mean relative refractory period
ANparams.cohc= 1.0;    % normal ohc function
ANparams.cihc= 1.0;    % normal ihc function
ANparams.species= 1;    % 1 for cat (2 for human with Shera et al. tuning; 3 for human with Glasberg & Moore tuning)
ANparams.noiseType= 0;  % 1 for variable fGn; 0 for fixed (frozen) fGn
ANparams.implnt= 0;     % "0" for approximate or "1" for actual implementation of the power-law functions in the Synapse
ANparams.dt= 1/StimParams.fs; %  time step

% Analysis params
anl.trifiltWidth= 5;
anl.onsetIgnore= 10e-3;
anl.fmSpread= 5;

total_nreps= 25;
useSpikes= 1;
useFilter= 1;
useFracPower= 1;

if useSpikes
    if useFilter
        postfix= 'LPed-Spikes';
    else
        postfix= 'Spikes';
    end
end


all_CF_Hz= logspace(log10(250), log10(8e3), 24); %[250 500 750 1.25e3 2.5e3 4e3 8e3];
example_CFs= [1e3 1.7e3 4e3];
ex_ind= dsearchn(all_CF_Hz(:), example_CFs(:));
all_CF_Hz(ex_ind)= example_CFs;

fsPlot= 40e3;

SCpeakAmp_org= nan(length(all_CF_Hz), 1);
SCpeakAmp_adj= nan(length(all_CF_Hz), 1);
PSTHavg_totPow= nan(length(all_CF_Hz), 1);
% PSTHeT_totPow= nan(length(all_CF_Hz), 1);
PSTHavg_fmPow_raw= nan(length(all_CF_Hz), 1);
PSTHeT_fmPow_raw= nan(length(all_CF_Hz), 1);
PSTHavg_fmPow_frac= nan(length(all_CF_Hz), 1);
PSTHeT_fmPow_frac= nan(length(all_CF_Hz), 1);
PSTHavg_RDow_raw= nan(length(all_CF_Hz), 1);
PSTHeT_RDow_raw= nan(length(all_CF_Hz), 1);

PSTHcT_LSB_SNR= nan(length(all_CF_Hz), 1);
PSTHcT_RSB_SNR= nan(length(all_CF_Hz), 1);
PSTHphiT_LSB_SNR= nan(length(all_CF_Hz), 1);
PSTHphiT_RSB_SNR= nan(length(all_CF_Hz), 1);

all_dBSPLused= nan(length(all_CF_Hz), 1);
all_drivenRates= nan(length(all_CF_Hz), 1);
VSpp_at_CF= nan(length(all_CF_Hz), 1);

ex_sT_PSD= cell(numel(all_CF_Hz), 1);
ex_eT_PSD= cell(numel(all_CF_Hz), 1);
ex_cT_PSD= cell(numel(all_CF_Hz), 1);
ex_phiT_PSD= cell(numel(all_CF_Hz), 1);
ex_freq_PSD= cell(numel(all_CF_Hz), 1);

rng('default');
% for freqVar= 10%:length(all_CF_Hz)
parfor freqVar= 1:length(all_CF_Hz)
    %% Create SAM tone and play it at 30 dB re Threshold
    curCF_Hz= all_CF_Hz(freqVar);
    
    
    [curSAM_pos, ~]= helper.create_SAM(curCF_Hz, StimParams.fm, StimParams.fs, StimParams.modDepth, StimParams.dur, [], StimParams.phi_m);
    
    %     thresh_dBSPL= get_thresh_curCF(curSAM_pos(1:round(StimParams.rlf_dur*StimParams.fs)), curCF_Hz, ANparams);
    %     curSAM_pos= gen_rescale(curSAM_pos, thresh_dBSPL + StimParams.dBreThresh);
    
    oa_dBSPL= helper.get_dBSPL_from_rlf(curSAM_pos(1:round(StimParams.rlf_dur*StimParams.fs)), curCF_Hz, ANparams, StimParams.DrivenRateTarget*StimParams.rlf_dur);
    %     oa_dBSPL= 30;

    curSAM_pos= helper.gen_rescale(curSAM_pos, oa_dBSPL);
    curSAM_neg= -curSAM_pos;
    
    %% Get spike trains
    nrep= total_nreps;
    tStart= 10e-3;
    tEnd= StimParams.dur;
    
    vihc_pos = ANModelBEZ2018.model_IHC_BEZ2018(curSAM_pos(:)' , curCF_Hz, nrep, ANparams.dt, 1.1*StimParams.dur, ANparams.cohc, ANparams.cihc,ANparams.species);
    [curSpikeTrainsPos, ~, ~, ~, ~, ~] = ANModelBEZ2018.model_Synapse_BEZ2018_trials(...
        vihc_pos, curCF_Hz, nrep, ANparams.dt, ANparams.noiseType, ANparams.implnt, ANparams.spont, ANparams.tabs, ANparams.trel);
    curSpikeTrainsPos= reshape(curSpikeTrainsPos, numel(curSpikeTrainsPos)/nrep, nrep);
    curSpikeTrainsPos= num2cell(curSpikeTrainsPos, 1);   
    curSpikeTrainsPos= cellfun(@(x,y) find(x)*y, curSpikeTrainsPos, repmat({ANparams.dt}, size(curSpikeTrainsPos)), 'UniformOutput', false);
    curSpikeTrainsPos= cellfun(@(x,t1,t2) x(x>t1 & x<t2), curSpikeTrainsPos, repmat({tStart}, size(curSpikeTrainsPos)), repmat({tEnd}, size(curSpikeTrainsPos)), 'UniformOutput', false);
    
    vihc_neg = ANModelBEZ2018.model_IHC_BEZ2018(curSAM_neg(:)' , curCF_Hz, nrep, ANparams.dt, 1.1*StimParams.dur, ANparams.cohc, ANparams.cihc,ANparams.species);
    [curSpikeTrainsNeg, ~, ~, ~, ~, ~] = ANModelBEZ2018.model_Synapse_BEZ2018_trials(...
        vihc_neg, curCF_Hz, nrep, ANparams.dt, ANparams.noiseType, ANparams.implnt, ANparams.spont, ANparams.tabs, ANparams.trel);
    curSpikeTrainsNeg= reshape(curSpikeTrainsNeg, numel(curSpikeTrainsNeg)/nrep, nrep);
    curSpikeTrainsNeg= num2cell(curSpikeTrainsNeg, 1);
    curSpikeTrainsNeg= cellfun(@(x,y) find(x)*y, curSpikeTrainsNeg, repmat({ANparams.dt}, size(curSpikeTrainsNeg)), 'UniformOutput', false);
    curSpikeTrainsNeg= cellfun(@(x,t1,t2) x(x>t1 & x<t2), curSpikeTrainsNeg, repmat({tStart}, size(curSpikeTrainsNeg)), repmat({tEnd}, size(curSpikeTrainsNeg)), 'UniformOutput', false);

    
    %% Compute SUMCOR peak
    [NSACpos,~,~,~] = helper.SACfull_m(curSpikeTrainsPos, 1/fsPlot, StimParams.dur);
    [NSACneg,~,~,~] = helper.SACfull_m(curSpikeTrainsNeg, 1/fsPlot, StimParams.dur);
    NSAC= helper.trifilt((NSACpos + NSACneg)/2, anl.trifiltWidth);
    [NSCC, ~,AVGrates,~] = helper.SCCfull_m({curSpikeTrainsPos, curSpikeTrainsNeg}, 1/fsPlot, StimParams.dur);
    NSCC= (NSCC + fliplr(NSCC))/2;
    NSCC= helper.trifilt(NSCC, anl.trifiltWidth);
    SUMCORraw= (NSAC + NSCC)/2;
    window_correction1= linspace(0, 1, (numel(SUMCORraw)+1)/2);
    window_correction= [window_correction1 fliplr(window_correction1(1:end-1))];
    SUMCOR= SUMCORraw - window_correction;
    SUMCOR_DC = mean(SUMCOR);
    DFT_SC= fft(SUMCOR - SUMCOR_DC);
    Freq_SC= linspace(0, fsPlot, numel(DFT_SC));

    SCzero_FreqMask= (Freq_SC > (StimParams.fm - anl.fmSpread) & Freq_SC < (StimParams.fm + anl.fmSpread)) ...
        | (Freq_SC > (fsPlot - StimParams.fm - anl.fmSpread) & Freq_SC < (fsPlot - StimParams.fm + anl.fmSpread)) ...
        | (Freq_SC > (2*StimParams.fm - anl.fmSpread) & Freq_SC < (2*StimParams.fm + anl.fmSpread)) ...
        | (Freq_SC > (fsPlot - 2*StimParams.fm - anl.fmSpread) & Freq_SC < (fsPlot - 2*StimParams.fm + anl.fmSpread)) ...
        | (Freq_SC > (3*StimParams.fm - anl.fmSpread) & Freq_SC < (3*StimParams.fm + anl.fmSpread)) ... 
        | (Freq_SC > (fsPlot - 3*StimParams.fm - anl.fmSpread) & Freq_SC < (fsPlot - 3*StimParams.fm + anl.fmSpread));

    SUMCORadj= real(ifft(DFT_SC .* SCzero_FreqMask)) + SUMCOR_DC;
    
    SCpeakAmp_org(freqVar)= max(SUMCOR);
    SCpeakAmp_adj(freqVar)= max(SUMCORadj);
    
    %%
    binEdges= 0:1/fsPlot:StimParams.dur;

    uRatePos= histcounts(cell2mat(curSpikeTrainsPos'), binEdges);
    uRateNeg= histcounts(cell2mat(curSpikeTrainsNeg'), binEdges);
    uRateAvg= (uRatePos + uRateNeg) / 2;
    uRateComp= (uRatePos - uRateNeg) / 2;
    
    if useFilter
        filtObj= helper.get_filter_fdesign('bp', [max(.1, curCF_Hz-5*StimParams.fm) curCF_Hz+5*StimParams.fm], fsPlot, 2); % second order for now
        uRateComp= filter(filtObj, uRateComp);
    end
    uRateHilbEnv= abs(hilbert(uRateComp))/sqrt(2);

    if useSpikes
        [Pxx_avg_dB, ~]= helper.plot_dpss_psd(uRateAvg, fsPlot, 'plot', false);
        hold on;
        [Pxx_eT_dB, FreqPSTH]= helper.plot_dpss_psd(uRateHilbEnv, fsPlot, 'plot', false);
        Pxx_cT_dB= helper.plot_dpss_psd(uRateComp, fsPlot, 'plot', false);
        Pxx_phiT_dB= helper.plot_dpss_psd( sqrt(2) * rms(uRateComp) * (uRateComp./abs(hilbert(uRateComp))), fsPlot, 'plot', false);
    end
    
    
    % Power near Fm
    %     warning('Change FM band');
    fm_mask_FreqPSTH= FreqPSTH > (StimParams.fm - anl.fmSpread) & FreqPSTH < (StimParams.fm + anl.fmSpread) ...
            | FreqPSTH > (2*StimParams.fm - anl.fmSpread) & FreqPSTH < (2*StimParams.fm + anl.fmSpread) ...
            | FreqPSTH > (3*StimParams.fm - anl.fmSpread) & FreqPSTH < (3*StimParams.fm + anl.fmSpread);
     PSTHavg_totPow(freqVar)= pow2db(rms(uRateAvg).^2);
    PSTHavg_fmPow_raw(freqVar)= pow2db(sum(db2pow(Pxx_avg_dB(fm_mask_FreqPSTH))));
    PSTHeT_fmPow_raw(freqVar)= pow2db(sum(db2pow(Pxx_eT_dB(fm_mask_FreqPSTH))));
    if useFracPower
        PSTHavg_fmPow_frac(freqVar)= db(sum(db2pow(Pxx_avg_dB(fm_mask_FreqPSTH))) / sum(db2pow(Pxx_avg_dB)));
        PSTHeT_fmPow_frac(freqVar)= db(sum(db2pow(Pxx_eT_dB(fm_mask_FreqPSTH))) / sum(db2pow(Pxx_eT_dB)));
    else
        PSTHavg_fmPow_frac(freqVar)= sqrt(sum(db2pow(Pxx_avg_dB(fm_mask_FreqPSTH)))) / mean(cellfun(@(x) mean(x), AVGrates));
        PSTHeT_fmPow_frac(freqVar)= sqrt(sum(db2pow(Pxx_eT_dB(fm_mask_FreqPSTH)))) / mean(cellfun(@(x) mean(x), AVGrates));
    end
    
    % Rectifier Distortion
    fcpmfm_mask_FreqPSTH= (FreqPSTH > (2*curCF_Hz - StimParams.fm - anl.fmSpread) & FreqPSTH < (2*curCF_Hz - StimParams.fm + anl.fmSpread)) ...
        | (FreqPSTH > (2*curCF_Hz - anl.fmSpread) & FreqPSTH < (2*curCF_Hz + anl.fmSpread)) ...
        | (FreqPSTH > (2*curCF_Hz + StimParams.fm - anl.fmSpread) & FreqPSTH < (2*curCF_Hz + StimParams.fm + anl.fmSpread));
    PSTHavg_RDow_raw(freqVar)= pow2db(sum(db2pow(Pxx_avg_dB(fcpmfm_mask_FreqPSTH))));
    PSTHeT_RDow_raw(freqVar)= pow2db(sum(db2pow(Pxx_eT_dB(fcpmfm_mask_FreqPSTH))));
    
    % SNR@Sidebands for TFS
    fc_m_fm_mask_FreqPSTH= (FreqPSTH > (curCF_Hz - StimParams.fm - anl.fmSpread) & FreqPSTH < (curCF_Hz - StimParams.fm + anl.fmSpread));
    fc_mask_FreqPSTH= (FreqPSTH > (curCF_Hz - anl.fmSpread) & FreqPSTH < (curCF_Hz + anl.fmSpread));
    fc_p_fm_mask_FreqPSTH= (FreqPSTH > (curCF_Hz + StimParams.fm - anl.fmSpread) & FreqPSTH < (curCF_Hz + StimParams.fm + anl.fmSpread));
    
    PSTHcT_LSB_SNR(freqVar)= pow2db(sum(db2pow(Pxx_cT_dB(fc_mask_FreqPSTH)))) - pow2db(sum(db2pow(Pxx_cT_dB(fc_m_fm_mask_FreqPSTH))));
    PSTHphiT_LSB_SNR(freqVar)= pow2db(sum(db2pow(Pxx_phiT_dB(fc_mask_FreqPSTH)))) - pow2db(sum(db2pow(Pxx_phiT_dB(fc_m_fm_mask_FreqPSTH))));
    PSTHcT_RSB_SNR(freqVar)= pow2db(sum(db2pow(Pxx_cT_dB(fc_mask_FreqPSTH)))) - pow2db(sum(db2pow(Pxx_cT_dB(fc_p_fm_mask_FreqPSTH))));
    PSTHphiT_RSB_SNR(freqVar)= pow2db(sum(db2pow(Pxx_phiT_dB(fc_mask_FreqPSTH)))) - pow2db(sum(db2pow(Pxx_phiT_dB(fc_p_fm_mask_FreqPSTH))));
    
    % Example PSDs
    all_dBSPLused(freqVar)= helper.calc_dbspl(curSAM_pos);
    all_drivenRates(freqVar)= mean(cellfun(@(x) mean(x), AVGrates));
    valid_CF_example= (ex_ind==freqVar);
    if any(valid_CF_example)
        ex_sT_PSD{freqVar}= Pxx_avg_dB;
        ex_eT_PSD{freqVar}= Pxx_eT_dB;
        ex_cT_PSD{freqVar}= Pxx_cT_dB;
        ex_phiT_PSD{freqVar}= Pxx_phiT_dB;
        ex_freq_PSD{freqVar}= FreqPSTH;
    else
        ex_sT_PSD{freqVar}= nan;
        ex_eT_PSD{freqVar}= nan;
        ex_cT_PSD{freqVar}= nan;
        ex_phiT_PSD{freqVar}= nan;
        ex_freq_PSD{freqVar}= nan;
    end
    
    % Get VSpp for tone@CF
    nrepTone= 30;
    tone_at_CF= helper.gen_rescale(helper.create_sinusoid(curCF_Hz, StimParams.fs, StimParams.dur), all_dBSPLused(freqVar));
    vihc_tone= ANModelBEZ2018.model_IHC_BEZ2018(tone_at_CF(:)' , curCF_Hz, 1, ANparams.dt, 1.1*StimParams.dur, ANparams.cohc, ANparams.cihc,ANparams.species);
    tonePST= cell(nrepTone, 1);
    for repVar= 1:nrepTone
        [tempPST, ~, ~, ~, ~, ~] = ANModelBEZ2018.model_Synapse_BEZ2018(...
            vihc_tone, curCF_Hz, 1, ANparams.dt, ANparams.noiseType, ANparams.implnt, ANparams.spont, ANparams.tabs, ANparams.trel);
        tonePST{repVar}= find(tempPST)*ANparams.dt;
    end
    VSpp_at_CF(freqVar)= helper.compute_VSpp(tonePST, curCF_Hz, anl.onsetIgnore, StimParams.dur);
end

%% ENV
figSize_cm= [15 5 19.05 13];
figure_prop_name = {'PaperPositionMode','units','Position', 'Renderer'};
figure_prop_val =  { 'auto'            ,'centimeters', figSize_cm, 'painters'};  % [Xcorner Ycorner Xwidth Ywidth]
figure(FigHan.ENV);
clf;
set(gcf,figure_prop_name,figure_prop_val);


plt.lw= 1;
plt.lw2= 1.25;
plt.mrkSize2= 7;
plt.mrkSize1= 5;
plt.fSize= 9;
plt.xtick_psd_val= [10 100 1e3 5e3];
plt.xtick_psd_lab= cellfun(@(x) num2str(x), num2cell(plt.xtick_psd_val/1e3), 'uniformoutput', false);

plt.xtick_CF_val= [300 500 1.5e3 3e3 5e3 10e3];
plt.xtick_CF_lab= cellfun(@(x) num2str(x), num2cell(plt.xtick_CF_val/1e3), 'uniformoutput', false);

plt.VSytick= [0 .5 1.0];
plt.SCpeakdB_ytick= floor(min(db(SCpeakAmp_adj))):4:ceil(max(db(SCpeakAmp_org)));

if ~useSpikes
    mrkYval2= 50;
    mrkYval1= 15;
else
    mrkYval2= -30;
    mrkYval1= -65;
end


panelLetters= 'ABC';

fmFactor1= 1.15;
fm_box1_X= StimParams.fm*[1/fmFactor1 fmFactor1 fmFactor1 1/fmFactor1 1/fmFactor1];
fmFactor2= 1.1;
fm_box2_X= 2*StimParams.fm*[1/fmFactor2 fmFactor2 fmFactor2 1/fmFactor2 1/fmFactor2];
fm_box3_X= 3*StimParams.fm*[1/fmFactor2 fmFactor2 fmFactor2 1/fmFactor2 1/fmFactor2];
fm_box_Y= [mrkYval1 mrkYval1 mrkYval2 mrkYval2 mrkYval1];

env_sp_ax= nan(8,1);
env_ttl_txt_Han= nan(8,1);

for exVar= 1:length(example_CFs)
    curCFind= ex_ind(exVar);
    tempCF_Hz= all_CF_Hz(curCFind);
    FreqPSTH= ex_freq_PSD{curCFind};
    
    env_sp_ax(exVar)= subplot(4, 3, exVar);
    hold on;
    plot(FreqPSTH, ex_sT_PSD{curCFind}, 'linew', plt.lw);
    plot(FreqPSTH, ex_eT_PSD{curCFind}, 'linew', plt.lw);
    
    % f_m box
    plot(fm_box1_X, fm_box_Y, '-.', 'color', helper.get_color('gray'), 'linew', plt.lw);
    plot(fm_box2_X, fm_box_Y, '-.', 'color', helper.get_color('gray'), 'linew', plt.lw);
    plot(fm_box3_X, fm_box_Y, '-.', 'color', helper.get_color('gray'), 'linew', plt.lw);

    % RD boxes    
    cf_box1_x1= .85*(2*tempCF_Hz - StimParams.fm);
    cf_box1_x2= 1.25*(2*tempCF_Hz - StimParams.fm);
    plot([cf_box1_x1 cf_box1_x2 cf_box1_x2 cf_box1_x1 cf_box1_x1], fm_box_Y, '-.', 'color', helper.get_color('br'), 'linew', plt.lw);
    
    % Markers for CF and set axes properties
    plot(tempCF_Hz, 50, 'v', 'color', 'k', 'linew', plt.lw);
    set(gca, 'XScale', 'log');
    set(gca,'xtick', plt.xtick_psd_val, 'XTickLabel', plt.xtick_psd_lab,'TickLength', plt.tick_len1, 'units', 'normalized'); 
    env_ttl_txt_Han(exVar)= text(0, 1.15, ['\bf' panelLetters(exVar) '\rm' sprintf('. CF=%.1f kHz', all_CF_Hz(curCFind)/1e3)], 'Units', 'normalized');
end
linkaxes(env_sp_ax(1:length(example_CFs)));
xlim([.7*StimParams.fm 10.2e3]);
if useSpikes
    ylim([mrkYval1-5 -30]);
else
    ylim([mrkYval1-2.5 50]);
end

[~, icons]= legend('S(f)', 'E(f)', 'box', 'off', 'Location', 'northeast');
icons(3).XData= mean(icons(3).XData) + [0 +.25];
icons(5).XData= mean(icons(5).XData) + [0 +.25];

set(gca, 'XScale', 'log');

subplot(4, 3, 1);
ylHan(1)= ylabel('PSD (dB/Hz)', 'Units', 'normalized');

subplot(4, 3, 2);
xlabel('Modulation Frequency (kHz)');

env_sp_ax(4)= subplot(412);
plot(all_CF_Hz, PSTHavg_fmPow_raw, '-s', 'linew', plt.lw, 'markersize', plt.mrkSize2);
hold on;
plot(all_CF_Hz, PSTHeT_fmPow_raw, '-o', 'linew', plt.lw, 'markersize', plt.mrkSize2);
set(gca, 'XScale', 'log');
if ~useSpikes
    ylim([0 60]);
end
% ylabel({'F_m-related'; 'power (dB)'})
ylHan(2)= ylabel('P_{Fm}(dB)', 'Units', 'normalized');
set(gca,'xtick', plt.xtick_CF_val, 'XTickLabel', plt.xtick_CF_lab,'TickLength', plt.tick_len2, 'units', 'normalized');

yyaxis right;
lHanVS= plot(all_CF_Hz, VSpp_at_CF, '--', 'color', helper.get_color('k'), 'linew', plt.lw, 'markersize', plt.mrkSize2);
set(gca, 'YColor', helper.get_color('k'), 'ytick', plt.VSytick);

ylim([-.01 1.01]);
env_ttl_txt_Han(4)= text(0.01, 1.1, '\bfD\rm', 'Units', 'normalized');

[~, icons]= legend(lHanVS, 'VS_p_p', 'box', 'off', 'Location', 'southwest');
icons(2).XData= mean(icons(2).XData) + [-.1 +.25];

env_sp_ax(5)= subplot(413);
plot(all_CF_Hz, PSTHavg_RDow_raw, '-s', 'linew', plt.lw, 'markersize', plt.mrkSize2);
hold on;
plot(all_CF_Hz, PSTHeT_RDow_raw, '-o', 'linew', plt.lw, 'markersize', plt.mrkSize2);

ylHan(3)= ylabel('P_{RD} (dB)', 'Units', 'normalized');
set(gca, 'XScale', 'log');
set(gca,'xtick', plt.xtick_CF_val, 'XTickLabel', plt.xtick_CF_lab,'TickLength', plt.tick_len2, 'units', 'normalized');
if ~useSpikes
    ylim([0 60]);
else 
    ylim([-42 -23]);
end
env_ttl_txt_Han(5)= text(0.01, 1.1, '\bfE\rm', 'Units', 'normalized');
xlabel('CF (kHz)');

yyaxis right;
plot(all_CF_Hz, VSpp_at_CF, '--', 'color', helper.get_color('k'), 'linew', plt.lw, 'markersize', plt.mrkSize2);
set(gca, 'YColor', helper.get_color('k'), 'ytick', plt.VSytick);
vs_ylHan= ylabel('VS_p_p (tone at CF)');
vs_ylHan.Position(1:2)= [10e3 1.25];
ylim([-.01 .91]);


env_sp_ax(6)= subplot(4,3,10);
plot(all_CF_Hz, db(SCpeakAmp_org), '-s', 'linew', plt.lw, 'markersize', plt.mrkSize1, 'color', helper.get_color('g'));
hold on;
plot(all_CF_Hz, db(SCpeakAmp_adj), '-o', 'linew', plt.lw, 'markersize', plt.mrkSize1, 'color', helper.get_color('prp'));
[lgHan, icons]= legend('Raw', 'Adj', 'box', 'off', 'Location', 'southeast');
icons(3).XData= mean(icons(3).XData) + [0 +.25];
icons(4).XData= icons(4).XData + .125;
icons(5).XData= mean(icons(5).XData) + [0 +.25];
icons(6).XData= icons(6).XData + .125;
lgHan.Position(1:2)= [.22 .095];

ylHan(4)= ylabel('SC_{peak} (dB)', 'Units', 'normalized');
set(gca, 'XScale', 'log');
set(gca,'xtick', plt.xtick_CF_val, 'XTickLabel', plt.xtick_CF_lab,'TickLength', plt.tick_len1, 'units', 'normalized', 'ytick', plt.SCpeakdB_ytick);
xlabel('CF (kHz)');
env_ttl_txt_Han(6)= text(0.05, 1.1, '\bfF\rm', 'Units', 'normalized');
ylim([min(plt.SCpeakdB_ytick)-.01 max(plt.SCpeakdB_ytick)+.01]);

linkaxes(env_sp_ax(4:6), 'x');
xlim([min(all_CF_Hz)*.95 1.05*max(all_CF_Hz)]);

env_sp_ax(7)= subplot(4,3,11);
plot(PSTHavg_fmPow_raw, db(SCpeakAmp_org), 's', 'linew', plt.lw, 'markersize', plt.mrkSize1, 'color', helper.get_color('g'));
hold on
plot(PSTHavg_fmPow_raw, db(SCpeakAmp_adj), 'o', 'linew', plt.lw, 'markersize', plt.mrkSize1, 'color', helper.get_color('prp'));

xlabel('P_{Fm} in s(t) (dB)');

env_ttl_txt_Han(7)= text(0.05, 1.1, '\bfG\rm', 'Units', 'normalized');
xlim([min(PSTHavg_fmPow_raw)-.2 max(PSTHavg_fmPow_raw)+.2]);
set(gca, 'TickLength', plt.tick_len1, 'units', 'normalized', 'ytick', plt.SCpeakdB_ytick);
ylim([min(plt.SCpeakdB_ytick)-.01 max(plt.SCpeakdB_ytick)+.01]);


env_sp_ax(8)= subplot(4,3,12);
plot(db2pow(PSTHavg_RDow_raw), SCpeakAmp_org - SCpeakAmp_adj, 'd', 'linew', plt.lw, 'markersize', plt.mrkSize1, 'color', helper.get_color('gray'))
xlabel('P_{RD}');
ylabel('SC_{peak}:Raw-Adj');
env_ttl_txt_Han(8)= text(0.05, 1.1, '\bfH\rm', 'Units', 'normalized');
xlim([min(db2pow(PSTHavg_RDow_raw))*.95 max(db2pow(PSTHavg_RDow_raw))*1.05]);
set(gca, 'TickLength', plt.tick_len1, 'units', 'normalized');
ylim([-.001 max(ylim)])

set(findall(gcf,'-property','fontsize'),'fontsize', plt.fSize);
set(env_ttl_txt_Han, 'FontSize', 11);
set(findall(gcf,'-property','box'),'box', 'off');

ylHan(1).Position(1)= -.15;
ylHan(2).Position(1)= -.04;
ylHan(3).Position(1)= -.04;
ylHan(4).Position(1)= -.15;

% ENV 
Xcorner_X= .07;
Xwidth_ABC= .255;
Xwidth_FGH= .24;
Xshift_ABC= .045;
Xshift_FGH= .06;
Xshift_ratio_FGH= .8;
Ycorner_X= .085;
Ywidth_X= .15;
Yshift_X= .09;

% F 
set(env_sp_ax(6),'Position',[Xcorner_X Ycorner_X Xwidth_FGH Ywidth_X])
drawnow

% G 
set(env_sp_ax(7),'Position',[Xcorner_X+Xwidth_FGH+Xshift_ratio_FGH*Xshift_FGH Ycorner_X Xwidth_FGH Ywidth_X])
drawnow

% H
set(env_sp_ax(8),'Position',[Xcorner_X+2*Xwidth_FGH+2*(2-Xshift_ratio_FGH)*Xshift_FGH Ycorner_X Xwidth_FGH Ywidth_X])
drawnow

% E
set(env_sp_ax(5),'Position',[Xcorner_X Ycorner_X+Ywidth_X+Yshift_X 3*Xwidth_ABC+2*Xshift_ABC Ywidth_X])
drawnow

% D
set(env_sp_ax(4),'Position',[Xcorner_X Ycorner_X+2*Ywidth_X+2*Yshift_X 3*Xwidth_ABC+2*Xshift_ABC Ywidth_X])
drawnow

% A
set(env_sp_ax(1),'Position',[Xcorner_X Ycorner_X+3*Ywidth_X+3*Yshift_X Xwidth_ABC Ywidth_X])
drawnow

% B
set(env_sp_ax(2),'Position',[Xcorner_X+Xwidth_ABC+Xshift_ABC Ycorner_X+3*Ywidth_X+3*Yshift_X Xwidth_ABC Ywidth_X])
drawnow

% C
set(env_sp_ax(3),'Position',[Xcorner_X+2*Xwidth_ABC+2*Xshift_ABC Ycorner_X+3*Ywidth_X+3*Yshift_X Xwidth_ABC Ywidth_X])
drawnow


if saveFig && (useSpikes)
    saveas(gcf, [LatexDir 'Fig6'], 'epsc');
end

%% TFS
figSize_cm= [15 5 13.2 9];
figure_prop_name = {'PaperPositionMode','units','Position', 'Renderer'};
figure_prop_val =  { 'auto'            ,'centimeters', figSize_cm, 'painters'};  % [Xcorner Ycorner Xwidth Ywidth]
figure(FigHan.TFS);
clf;
set(gcf,figure_prop_name,figure_prop_val);

if ~useSpikes
    mrkYval2= 55;
else
    mrkYval2= -25;
end

figure(FigHan.TFS);
clf;
tfs_sp_ax= nan(5, 1);
tfs_ttl_txt_Han= nan(5,1);
for exVar= 1:length(example_CFs)
    curCFind= ex_ind(exVar);
    tempCF_Hz= all_CF_Hz(curCFind);
    FreqPSTH= ex_freq_PSD{curCFind};
    
    cur_xtick_psd_val= example_CFs(exVar)+StimParams.fm*(-6:3:6);
    cur_xtick_psd_lab= cellfun(@(x) num2str(x), num2cell(cur_xtick_psd_val), 'uniformoutput', false);
    
    
    tfs_sp_ax(exVar)= subplot(3, 3, exVar);
    hold on;
    plot(FreqPSTH, ex_cT_PSD{curCFind}, 'linew', plt.lw);
    plot(FreqPSTH, ex_phiT_PSD{curCFind}, 'LineStyle', '--', 'linew', plt.lw);
    plot(tempCF_Hz, mrkYval2, 'v', 'color', 'k', 'linew', plt.lw);
    plot(tempCF_Hz - StimParams.fm, mrkYval2, 'o', 'color', helper.get_color('prp'), 'linew', plt.lw);
    plot(tempCF_Hz + StimParams.fm, mrkYval2, 's', 'color', helper.get_color('prp'), 'linew', plt.lw);
    set(gca, 'XScale', 'log');
    tfs_ttl_txt_Han(exVar)= text(0.05, 1.15, ['\bf' panelLetters(exVar) '\rm' sprintf('. CF=%.1f kHz', all_CF_Hz(curCFind)/1e3)], 'Units', 'normalized');
    set(gca,'xtick', cur_xtick_psd_val, 'XTickLabel', cur_xtick_psd_lab);
    xlim(tempCF_Hz + [-5*StimParams.fm 5*StimParams.fm]);
    set(gca, 'TickLength', plt.tick_len1, 'units', 'normalized');
    if ~useSpikes
        set(gca, 'YTick', 10*(round((mrkYval2-40)/10):round((mrkYval2+2)/10)));
    end
end
linkaxes(tfs_sp_ax(1:3), 'y');
if ~useSpikes
    ylim([0 60]);
else 
    ylim([mrkYval2-40 mrkYval2+2]);
end

[~, icons]= legend('D(f)', '\Phi(f)', 'box', 'off', 'Location', 'northeast');
icons(3).XData= mean(icons(3).XData) + [-.1 +.3];
icons(4).XData= icons(4).XData + .125;
icons(5).XData= mean(icons(5).XData) + [-.1 +.3];
icons(6).XData= icons(6).XData + .125;
set(gca, 'XScale', 'log');

subplot(3, 3, 1);
ylabel('PSD (dB/Hz)');

subplot(3, 3, 2);
xlabel('Carrier Frequency (Hz)');

tfs_sp_ax(4)= subplot(312);
plot(all_CF_Hz, PSTHcT_LSB_SNR, '-s', 'linew', plt.lw, 'markersize', plt.mrkSize2);
hold on;
plot(all_CF_Hz, PSTHphiT_LSB_SNR, '-o', 'linew', plt.lw, 'markersize', plt.mrkSize2);
set(gca, 'XScale', 'log');
ylabel('P_{F_c}/P_{LSB} (dB)');
set(gca,'xtick', plt.xtick_CF_val, 'XTickLabel', plt.xtick_CF_lab);
tfs_ttl_txt_Han(4)= text(0.01, 1.1, '\bfD\rm', 'Units', 'normalized');
xlim([min(all_CF_Hz)*.95 1.05*max(all_CF_Hz)]);
set(findall(gca,'-property','TickLength'),'TickLength', plt.tick_len2, 'units', 'normalized');

tfs_sp_ax(5)= subplot(313);
plot(all_CF_Hz, PSTHcT_RSB_SNR, '-s', 'linew', plt.lw, 'markersize', plt.mrkSize2);
hold on;
plot(all_CF_Hz, PSTHphiT_RSB_SNR, '-o', 'linew', plt.lw, 'markersize', plt.mrkSize2);
set(gca, 'XScale', 'log');
ylabel('P_{F_c}/P_{USB} (dB)');
set(gca,'xtick', plt.xtick_CF_val, 'XTickLabel', plt.xtick_CF_lab);
legend('D(f)', '\Phi(f)', 'box', 'off', 'Location', 'northeast');
tfs_ttl_txt_Han(5)= text(0.01, 1.1, '\bfE\rm', 'Units', 'normalized');
xlabel('CF (kHz)');
xlim([min(all_CF_Hz)*.95 1.05*max(all_CF_Hz)]);

set(findall(gcf,'-property','fontsize'),'fontsize', plt.fSize);
set(findall(gcf,'-property','box'),'box', 'off');
set(findall(gca,'-property','TickLength'),'TickLength', plt.tick_len2, 'units', 'normalized');
set(tfs_ttl_txt_Han, 'FontSize', 11);

% TFS
Xcorner_X= .095;
Xwidth_ABC= .265;
Xshift_ABC= .05;
Ycorner_X= .12;
Ywidth_X= .205;
Yshift_X= .1;

% E
set(tfs_sp_ax(5),'Position',[Xcorner_X Ycorner_X 3*Xwidth_ABC+2*Xshift_ABC Ywidth_X])
drawnow

% D
set(tfs_sp_ax(4),'Position',[Xcorner_X Ycorner_X+Ywidth_X+Yshift_X 3*Xwidth_ABC+2*Xshift_ABC Ywidth_X])
drawnow

% A
set(tfs_sp_ax(1),'Position',[Xcorner_X Ycorner_X+2*Ywidth_X+2*Yshift_X Xwidth_ABC Ywidth_X])
drawnow

% B
set(tfs_sp_ax(2),'Position',[Xcorner_X+Xwidth_ABC+Xshift_ABC Ycorner_X+2*Ywidth_X+2*Yshift_X Xwidth_ABC Ywidth_X])
drawnow

% C
set(tfs_sp_ax(3),'Position',[Xcorner_X+2*Xwidth_ABC+2*Xshift_ABC Ycorner_X+2*Ywidth_X+2*Yshift_X Xwidth_ABC Ywidth_X])
drawnow


if saveFig && (useSpikes)
    saveas(gcf, [LatexDir 'Fig7'], 'epsc');
end
