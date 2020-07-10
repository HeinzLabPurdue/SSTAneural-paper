function Hd= get_filter_fdesign(filtType, freqWindow, fs, filtOrder, plotYes, verbose)

if ~exist('filtOrder', 'var')
    filtOrder= 4;
end
if ~exist('plotYes', 'var')
    plotYes =0;
end
if ~exist('verbose', 'var')
    verbose= 0;
end

switch filtType
    case {'band', 'bandpass', 'bp'}
        filtObject= fdesign.bandpass('N,F3dB1,F3dB2', filtOrder, freqWindow(1), freqWindow(2), fs);
        
    case {'low', 'lowpass', 'lp'}
        filtObject= fdesign.lowpass('N,F3dB', filtOrder, freqWindow, fs);
        
    case {'high', 'highpass', 'hp'}
        filtObject= fdesign.lowpass('N,F3dB', filtOrder, freqWindow, fs);
        
end

availabel_designMethods= designmethods(filtObject);
Hd= design(filtObject, availabel_designMethods{1});
if verbose
    fprintf('Using %s as design method \n', availabel_designMethods{1});
end


if plotYes
    freqz(Hd);
end