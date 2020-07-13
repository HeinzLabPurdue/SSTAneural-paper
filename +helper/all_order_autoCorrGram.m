function [NSAC,delay,AVGrate,TOTALspikes] = all_order_autoCorrGram(SpikeTrain1,DELAYbinwidth,Duration)
% File: SAChalf
% 21 Jun, 2017: SP
% Calls: SAChalf.c (mex)

% Computes Normalized Shuffled Auto-Coprrelogram (NSAC) from a set of Spike Trains and Duration
% Based on Louage et. al 2004


NUMspikeREPS=length(SpikeTrain1);
%%%%%%%%%%%% Compute AVGrate
Kmax=0;
for spikeREPind=1:NUMspikeREPS
    if length(SpikeTrain1{spikeREPind})>Kmax
        Kmax=length(SpikeTrain1{spikeREPind});
    end
end

%%%%%%%%%%%%% Setup BIG spike matrix for faster computation
SpikeMAT=NaN*ones(NUMspikeREPS,Kmax);
for REPindREF=1:NUMspikeREPS
    SpikeMAT(REPindREF,1:length(SpikeTrain1{REPindREF}))=SpikeTrain1{REPindREF};
end

NUMspikes=sum(~isnan(SpikeMAT(:,:))',1);  % Count number of real spikes in each line
TOTALspikes=sum(NUMspikes);
AVGrate=TOTALspikes/NUMspikeREPS/Duration;

%%%%%%%%%%%%% Compute Shuffled Auto-Correlogram
% SpikeMAT1 is setup in Matlab as: rows hold each spike train;
% C/MEX: indexing goes down 1st column 1st, so we need to pass the transpose of SpikeMAT1, to get easy indexing in MEXfile


[SAC,delay] = helper.AllOrder_autoCorrGram(SpikeMAT',NUMspikes,TOTALspikes,Duration, DELAYbinwidth);

% delay=[-fliplr(delay(2:end)) delay];
% SAC=[fliplr(intsMEX(2:end)) intsMEX];

NSAC=SAC; %/(NUMspikeREPS*(NUMspikeREPS-1)*AVGrate^2*DELAYbinwidth*Duration);  % From Louage et al (2004: J. Neurophysiol)

return;