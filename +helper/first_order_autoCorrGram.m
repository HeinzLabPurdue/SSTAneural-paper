function [foISIhist,delay] = first_order_autoCorrGram(SpikeTrain1,DELAYbinwidth,Duration)
% File: SAChalf
% 21 Jun, 2017: SP
% Calls: SAChalf.c (mex)

% Computes Normalized Shuffled Auto-Coprrelogram (NSAC) from a set of Spike Trains and Duration
% Based on Louage et. al 2004


NUMspikeREPS=length(SpikeTrain1);
histEdges= DELAYbinwidth/2:DELAYbinwidth:(Duration+DELAYbinwidth/2);
histEdges= [-fliplr(histEdges) histEdges];
delay= (histEdges(1:end-1)+histEdges(2:end))/2;
foISIhist= zeros(NUMspikeREPS, length(delay));

for repVar=1:NUMspikeREPS
    curSpikes= SpikeTrain1{repVar};
    first_order_ISI= diff(curSpikes);
    temp_foISI= histcounts(first_order_ISI, histEdges);
    foISIhist(repVar, :)= (temp_foISI + fliplr(temp_foISI))/2;
end

foISIhist= sum(foISIhist, 1);