function [NSCC,delay,AVGrates,TOTALspikes] = SCCfull_m(SpikeTrains,DELAYbinwidth,Duration)
% function [NSCC,delay,AVGrates,TOTALspikes] = SCCfull_m(SpikeTrains,DELAYbinwidth,Duration)
% File: SCCfull_m
% 22Jun2017: SP: Updated ShufCrossCorr
% 
% Calls: SCCfull.c (mex)
% Computes Normalized Shuffled Cross-Correlogram (NSCC) from a set of Spike Trains and Duration
% Based on Louage et. al 2004


%% Compute AVGrates
NUMspikeREPS=cell(1,2);
Kmax=cell(1,2);
NUMspikes=cell(1,2);
TOTALspikes=cell(1,2);
AVGrates=cell(1,2);
SpikeMAT=cell(1,2);

%%
for UNITind=1:2
    NUMspikeREPS{UNITind}=length(SpikeTrains{UNITind});
    NUMspikes{UNITind}=cellfun(@(x) numel(x),SpikeTrains{UNITind})';
    Kmax{UNITind}=max(NUMspikes{UNITind});
    TOTALspikes{UNITind}=sum(NUMspikes{UNITind});
    AVGrates{UNITind}=TOTALspikes{UNITind}/NUMspikeREPS{UNITind}/Duration;
    SpikeMAT{UNITind}=NaN*ones(NUMspikeREPS{UNITind},Kmax{UNITind});
    for REPindREF=1:NUMspikeREPS{UNITind}
        SpikeMAT{UNITind}(REPindREF,1:length(SpikeTrains{UNITind}{REPindREF}))=SpikeTrains{UNITind}{REPindREF};
    end    
end

[SCC,delay] = helper.SCCfull(SpikeMAT{1}',NUMspikes{1},TOTALspikes{1},SpikeMAT{2}',NUMspikes{2},TOTALspikes{2}, Duration, DELAYbinwidth);

NSCC=SCC/(NUMspikeREPS{1}*NUMspikeREPS{2}*AVGrates{1}*AVGrates{2}*DELAYbinwidth*Duration);

return;