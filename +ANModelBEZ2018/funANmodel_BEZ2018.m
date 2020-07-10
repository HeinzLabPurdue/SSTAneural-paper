function spikes_trains= funANmodel_BEZ2018(stim, fsStim, CF_Hz, spont, nrep , doPlot)
stim= stim(:)';

% model parameters
% CF_Hz    = 1e3;   % CF in Hz;

if ~exist('spont', 'var')
    spont = 100;   % spontaneous firing rate
end

if ~exist('doPlot', 'var')
    doPlot= false;
end

if ~exist('nrep', 'var')
    nrep = 1;               % number of stimulus repetitions (e.g., 50);
end


tabs   = 0.6e-3; % Absolute refractory period
trel   = 0.6e-3; % Baseline mean relative refractory period
cohc  = 1.0;    % normal ohc function
cihc  = 1.0;    % normal ihc function
species = 1;    % 1 for cat (2 for human with Shera et al. tuning; 3 for human with Glasberg & Moore tuning)
noiseType = 1;  % 1 for variable fGn; 0 for fixed (frozen) fGn
implnt = 0;     % "0" for approximate or "1" for actual implementation of the power-law functions in the Synapse

% stimulus parameters
% stimdb = 60; % stimulus intensity in dB SPL
% F0 = CF_Hz;     % stimulus frequency in Hz
Fs_AN = fsStim;  % sampling rate in Hz (must be 100, 200 or 500 kHz)
T  = numel(stim)/fsStim;  % stimulus duration in seconds
% rt = 2.5e-3; % rise/fall time in seconds
% ondelay = 10e-3;

% PSTH parameters

% t = 0:1/Fs:T-1/Fs; % time vector
% mxpts = length(t);
% irpts = rt*Fs;
% onbin = round(ondelay*Fs);

% pin = zeros(1,onbin+mxpts);
% pin(onbin+1:onbin+mxpts) = sqrt(2)*20e-6*10^(stimdb/20)*sin(2*pi*F0*t); % unramped stimulus
% pin(onbin+1:onbin+irpts)= pin(onbin+1:onbin+irpts).*(0:(irpts-1))/irpts;
% pin(onbin+(mxpts-irpts):onbin+mxpts)=pin(onbin+(mxpts-irpts):onbin+mxpts).*(irpts:-1:0)/irpts;

dt=1/Fs_AN; %  time step
uRate_duration= 1.05*T;
vihc = ANModelBEZ2018.model_IHC_BEZ2018(stim,CF_Hz,1,dt,uRate_duration,cohc,cihc,species);
vihc_nRep= repmat(vihc, 1, nrep);
[psth_nRep,meanrate, varrate, synout, trd_vector,trel_vector] = ANModelBEZ2018.model_Synapse_BEZ2018(vihc_nRep,CF_Hz,1,dt,noiseType,implnt,spont,tabs,trel);
psth= reshape(psth_nRep, numel(vihc), nrep);
psth_cell= num2cell(psth, 1)';
Fs_AN_cell= repmat({Fs_AN}, nrep, 1);
uRate_dur_cell= repmat({uRate_duration}, nrep, 1);
spikes_trains= cellfun(@(x,y,z) mod(find(x)/y, z), psth_cell, Fs_AN_cell, uRate_dur_cell, 'UniformOutput', false);
psth= nanmean(psth, 2)';

if doPlot
    psthbinwidth = 1e-4; % binwidth in seconds;
    psthbins = round(psthbinwidth*Fs_AN);  % number of psth bins per psth bin
    Psth = sum(reshape(psth,psthbins,length(psth)/psthbins)); %
    simtime = length(psth)/Fs_AN;
    tvect = 0:psthbinwidth:simtime-psthbinwidth;
    
    tt= 0:1/Fs_AN:(length(psth)-1)/Fs_AN;
    
    px = zeros(size(psth));
    px(1:length(stim)) = stim;
    
    Sout = mean(reshape(synout,length(synout)/nrep,nrep),2); %
    T_rd = mean(reshape(trd_vector,length(trd_vector)/nrep,nrep),2); %
    T_rel = mean(reshape(trel_vector,length(trel_vector)/nrep,nrep),2); %
    
    
    figure(133);
    subplot(3,1,1)
    plot(tt*1e3,px)
    ylabel('Pressure (Pa)')
    xlabel('Time (ms)')
    xlim(ceil(tt([1 end])*1e3))
    title('Acoustic Stimulus')
    subplot(3,1,2)
    plot(tt*1e3,vihc(1:length(tt))*1e3)
    ylabel('V_{ihc} (mV)')
    xlabel('Time (ms)')
    title('IHC relative membrane potential')
    xlim(ceil(tt([1 end])*1e3))
    subplot(3,1,3)
    bar(tvect*1e3, Psth/nrep/psthbinwidth,'histc') % Plot of estimated mean spike rate
    ylabel('Firing Rate (/s)')
    xlabel('Time (ms)')
    xlim(ceil(tt([1 end])*1e3))
    title('PSTH')
    
%     figure
%     subplot(3,1,1)
%     plot(tt*1e3,Sout)
%     ylabel('S_{out} (/s)')
%     xlabel('Time (ms)')
%     xlim(ceil(tt([1 end])*1e3))
%     title('Mean Synaptic Release Rate')
%     subplot(3,1,2)
%     plot(tt*1e3,meanrate)
%     ylabel('Mean Rate (/s)')
%     xlabel('Time (ms)')
%     title('Mean of Spike Rate')
%     xlim(ceil(tt([1 end])*1e3))
%     subplot(3,1,3)
%     plot(tt*1e3,varrate)
%     ylabel('Var Rate (/s)')
%     xlabel('Time (ms)')
%     xlim(ceil(tt([1 end])*1e3))
%     title('Variance in Spike Rate')
%     
%     figure;
%     subplot(3,1,1)
%     plot(tt*1e3,Sout)
%     ylabel('S_{out} (/s)')
%     xlabel('Time (ms)')
%     xlim(ceil(tt([1 end])*1e3))
%     title('Mean Synaptic Release Rate')
%     subplot(3,1,2)
%     plot(tt*1e3,T_rd*1e3)
%     ylabel('\tau_{rd} (ms)')
%     xlabel('Time (ms)')
%     title('Mean Synaptic Redocking Time')
%     xlim(ceil(tt([1 end])*1e3))
%     subplot(3,1,3)
%     plot(tt*1e3,T_rel*1e3)
%     ylabel('t_{rel} (ms)')
%     xlabel('Time (ms)')
%     xlim(ceil(tt([1 end])*1e3))
%     title('Mean Relative Refractory Period')
end