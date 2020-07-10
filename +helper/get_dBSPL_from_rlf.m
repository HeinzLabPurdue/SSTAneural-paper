function oa_dBSPL= get_dBSPL_from_rlf(stim, CF_Hz, ANparams, target_rate)

all_dBSPL= -20:3:100;
all_nSpikes= nan(size(all_dBSPL));

anl.nrep= 1;
anl.dur= length(stim)*ANparams.dt;
anl.onsetIgnore= 10e-3;

for splVar= 1:length(all_dBSPL)
    stim= helper.gen_rescale(stim, all_dBSPL(splVar));
    vihc_pos = ANModelBEZ2018.model_IHC_BEZ2018(stim(:)' , CF_Hz, anl.nrep, ANparams.dt, 1.1*anl.dur, ANparams.cohc, ANparams.cihc,ANparams.species);
    [psth_pos, ~, ~, ~, ~, ~] = ANModelBEZ2018.model_Synapse_BEZ2018(...
        vihc_pos, CF_Hz, anl.nrep, ANparams.dt, ANparams.noiseType, ANparams.implnt, ANparams.spont, ANparams.tabs, ANparams.trel);
    spike_times= find(psth_pos)*ANparams.dt;
    all_nSpikes(splVar)= sum(spike_times>anl.onsetIgnore & spike_times<anl.dur);
end

plotRLF= 0;
if plotRLF
    figure(2347);
    hold on;
    % plot(all_dBSPL, all_nSpikes, '-*');
end

plotFittedRLV= plotRLF;

[RLVparams,~,~]= helper.fitRLfun(all_dBSPL,all_nSpikes,plotRLF,plotFittedRLV);

oa_dBSPL= all_dBSPL(dsearchn(RLVparams.rates_est(:), target_rate));