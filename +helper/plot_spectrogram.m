function [timePlot, freqPlot, dB_powPlot]= plot_spectrogram(sig, fs, tWindow, fracOVlap, nfft, useDefaultPlot, tStart)
if ~exist('tStart', 'var')
    tStart= 0;
end

if ~exist('useDefaultPlot', 'var')
    useDefaultPlot= 1;
end
if ~exist('tWindow', 'var')
    tWindow= 64e-3;
elseif isempty(tWindow)
    tWindow= 64e-3;
end
if ~exist('fracOVlap', 'var')
    fracOVlap= .99;
elseif isempty(fracOVlap)
    fracOVlap= .99;
end

if ~exist('nfft', 'var')
    nfft= 2^(nextpow2(tWindow*fs)+1);
elseif isempty(nfft)
    nfft= 2^(nextpow2(tWindow*fs)+1);
end
nWindow= round(tWindow*fs);
if useDefaultPlot
    spectrogram(sig, blackman(nWindow), round(nWindow*fracOVlap), nfft, 'yaxis', fs);
    timePlot= nan;
    freqPlot= nan;
    dB_powPlot= nan;
    
else
    [S,F,T]= spectrogram(sig, blackman(nWindow), round(nWindow*fracOVlap), nfft, 'yaxis', fs);
    timePlot= (tStart+T)*1e3;
    freqPlot= F/1e3;
    dB_powPlot= pow2db(abs(S));
    surf(timePlot, freqPlot, dB_powPlot,'EdgeColor','none');
    view(2);
    grid off;
end