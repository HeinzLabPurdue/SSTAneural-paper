function [splVals, timeVals]= gen_get_spl_vals(sig, fs, tWindow, fracOverlap)

if ~exist('fracOVlap', 'var')
    fracOverlap= 0;
end

t= (1:length(sig))/fs;
tOVlap= tWindow*fracOverlap;
tSlide= tWindow-tOVlap;
pRef= 20e-6;

nSegs= 1 + round((t(end)-tWindow)/tSlide);

timeVals= nan(nSegs, 1);
splVals= nan(nSegs, 1);

for segVar= 1:nSegs
    seg_t_start=  (segVar-1)*tSlide;
    seg_ind_start= max(1, round(seg_t_start*fs)+1);
    seg_t_end= seg_t_start + tWindow;
    seg_ind_end= min(length(sig), round(seg_t_end*fs));
    
   timeVals(segVar)= (seg_t_start+seg_t_end)/2; 
   temp_sig= sig(seg_ind_start:seg_ind_end);
   rmsVal= rms(temp_sig-mean(temp_sig));
   splVals(segVar)= 20*log10(rmsVal/pRef);
end