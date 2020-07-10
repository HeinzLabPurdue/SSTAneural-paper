function vecOut= gen_rescale(vecIn, newSPL, verbose)

if ~exist('verbose', 'var')
    verbose=0;
end

if ~ismember(size(vecIn,2), [1,2])
    error('signal should be a one or two column matrix');
end

pRef= 20e-6; % for re. dB SPl
vecOut= nan(size(vecIn));
for chanVar= 1:size(vecIn, 2)
    if any(vecIn(:,chanVar))
        oldSPL= 20*log10(rms(vecIn(:,chanVar))/pRef);
        gainVal= 10^( (newSPL-oldSPL)/20 );
    else
        gainVal =0;
    end
    vecOut(:,chanVar)= vecIn(:,chanVar)*gainVal;
end

if verbose
    fprintf('Signal RMS= %.1f (Desired %.1f) \n', 20*log10(rms(vecOut)/pRef), newSPL);
end
