% Tested!!
% Filter a signal along a spectrotemporal trajectory. numel(inSig) has
% to be equal to numel(freqTrajectory).
% function fracSignal= get_trajectory_signal(inSig, fs, freqTrajectory, freq_half_window_co_Hz)
function [fracSignal, filtSignal]= get_trajectory_signal(inSig, fs, freqTrajectory, filtParams)

%% Make everything column vectors
inSig= inSig(:);
freqTrajectory= freqTrajectory(:);
freqTrajectory(isnan(freqTrajectory))= 0;
siglen= length(inSig);
stim_dur= siglen/fs;


if ~exist('filtParams', 'var') % doesn't exist
    filtParams= 2/stim_dur;
    d_lp = designfilt('lowpassiir','FilterOrder', 3, ...
        'HalfPowerFrequency', filtParams/(fs/2), 'DesignMethod','butter');
elseif isempty(filtParams) % exists but is empty
    filtParams= 2/stim_dur;
    d_lp = designfilt('lowpassiir','FilterOrder', 3, ...
        'HalfPowerFrequency', filtParams/(fs/2), 'DesignMethod','butter');
elseif isnumeric(filtParams) % exists and is numeric 
    d_lp = designfilt('lowpassiir','FilterOrder', 3, ...
        'HalfPowerFrequency', filtParams/(fs/2), 'DesignMethod','butter');
else % is a filter 
    d_lp= filtParams;
end


%%
phi_trajectory= -cumtrapz(freqTrajectory)/fs;
sig_demod_empirical= inSig .* exp(2*pi*1j*phi_trajectory); % sqrt_2 times to conserve power after Hilbert transform

filtSignal= sqrt(2)*filtfilt(d_lp,sig_demod_empirical); % multiply by sqrt(2) because have only considered +ve frequencies,
% have to consider for -ve frequencies too => power should be twice => amplitdue should be sqrt(2) times

filtSignal= abs(filtSignal); 
fracSignal= abs(filtSignal) / rms(inSig)^2; % normalize by total power