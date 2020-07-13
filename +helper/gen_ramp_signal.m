% function sigOUT= gen_ramp_signal(sigIN, fs, durONramp, durOFFramp, verbose)
function sigOUT= gen_ramp_signal(sigIN, fs, durONramp, durOFFramp, verbose)

if ~exist('verbose', 'var')
    verbose=false;
end
if nargin <3
    error('Not enough input parameters');
end
if ~exist('durOFFramp', 'var')
    durOFFramp= durONramp;
end

% create ON ramp
nSampON= round(fs*durONramp);
NW=1;


[ramp_tapers_on_all, ramp_weights_on_all]= dpss(2*nSampON, NW);
on_taper= ramp_tapers_on_all(1:nSampON,1);
if durOFFramp == durONramp
    off_taper=ramp_tapers_on_all((nSampON+1):end,1);
    ramp_weights_off_all= ramp_weights_on_all;
else
    nSampOFF= round(fs*durOFFramp);
    [ramp_tapers_off, ramp_weights_off_all]= dpss(2*nSampOFF, NW);
    off_taper= ramp_tapers_off((nSampOFF+1):end,1);
end

on_taper= on_taper/max(on_taper);
off_taper=off_taper/max(off_taper);

if verbose
    fprintf('Weights of On =%.3f and OFF=%.3f\n', ramp_weights_on_all(1), ramp_weights_off_all(1));
end

nOnes= length(sigIN)-length(on_taper)-length(off_taper);
ramp_vector= [on_taper; ones(nOnes, 1); off_taper];

if size(sigIN, 1)==1
    sigOUT= sigIN .*  ramp_vector';
elseif size(sigIN,2)==1
    sigOUT= sigIN .*  ramp_vector;
else
    error('Haven''t thought what to do with you!');
end