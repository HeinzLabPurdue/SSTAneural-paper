% Need to test more. 
% Compute power in a signal along a frequency trajectory. numel(inSig) has
% to be equal to numel(freqTrajectory). Options to use both multitaper PSD
% (pmtm) or DFT.
function [outPower, totPower, Pxx_dB, Freq_xx]= get_freq_trajectory_power(inSig, fs, freqTrajectory, plotPSD, freq_half_window_co_Hz, stim_freq_shift, nw, use_DFT0_DPSS1, freqMin)

if ~exist('freqMin', 'var')
    freqMin= -inf; % minimum frequency to compute total power
end

if ~exist('use_DFT0_DPSS1', 'var')
    use_DFT0_DPSS1= 1;
end
if ~exist('plotPSD', 'var')
    plotPSD= false;
end

%% Make everything column vectors
inSig= inSig(:);
freqTrajectory= freqTrajectory(:);

%%
siglen= length(inSig);
stim_dur= siglen/fs;
if ~exist('freq_half_window_co_Hz', 'var')
    freq_half_window_co_Hz= 2/stim_dur;
elseif isempty(freq_half_window_co_Hz)
    freq_half_window_co_Hz= 2/stim_dur;
end

stim_t= ((1:siglen)/fs)';
if ~exist('stim_freq_shift', 'var')
    stim_freq_shift= geomean(freqTrajectory);
end


phi_trajectory= -cumtrapz(freqTrajectory)/fs;
sig_demod_empirical= hilbert(inSig)/sqrt(2) .* exp(2*pi*1j*phi_trajectory); % sqrt_2 times to conserve power after Hilbert transform

sig_Plot= sig_demod_empirical.*exp(2*pi*1j*stim_freq_shift*stim_t); % Easier to plot PSD of FMed y to a different frequency than plotting DC.

if use_DFT0_DPSS1==0
    dft_sidedness= 1;
    [Adb, F_dft]= helper.plot_dft(sig_Plot, fs, 'sided', dft_sidedness, 'dc', true, 'plot', plotPSD==1);
    dft_validFreqInds= F_dft>(stim_freq_shift-freq_half_window_co_Hz) & F_dft<(stim_freq_shift+freq_half_window_co_Hz);
    dft_legalFreqInds= F_dft>freqMin;
    if plotPSD
        hold on;
        plot([stim_freq_shift-freq_half_window_co_Hz stim_freq_shift+freq_half_window_co_Hz], [max(Adb)+2 max(Adb)+2], 'color', helper.get_color('g'), 'LineWidth', 4);
    end
    
    if dft_sidedness==1
        dft_mask= 2*ones(size(Adb));
        dft_mask(1)= 1;
        nfft_true= numel(Adb)*2;
    elseif dft_sidedness==2
        dft_mask= ones(size(Adb));
        nfft_true= numel(Adb);
    end
      
    outPower= sum((db2mag(Adb(dft_validFreqInds))./dft_mask(dft_validFreqInds)).^2)*numel(sig_demod_empirical)/nfft_true;
    totPower= sum((db2mag(Adb(dft_legalFreqInds))./dft_mask(dft_legalFreqInds)).^2)*numel(sig_demod_empirical)/nfft_true;
    Pxx_dB= (db2mag(Adb(dft_validFreqInds))./dft_mask(dft_validFreqInds)).^2*numel(sig_demod_empirical)/nfft_true;
    Freq_xx= F_dft;
else
    
    if ~exist('nw', 'var')
        nw= 1.5;
    elseif isempty(nw)
        nw= 1.5;
    end
    [Pxx_dB, Freq_xx]= helper.plot_dpss_psd(sqrt(2)*real(sig_Plot), fs, 'nw', nw, 'title', '', 'norm', true, 'plot', plotPSD==1, 'yrange', 60);
    % See Table 8.1 Oppenheim and Schafer. DFT(Real(X)) = 0.5 (DFT(X) +
    % DFT(X*)). So to conserve power (after getting rid of negative
    % frequencies using Hilbert), multiply by sqrt(2).
    
    if plotPSD
        hold on;
        plot([stim_freq_shift-freq_half_window_co_Hz stim_freq_shift+freq_half_window_co_Hz], [max(Pxx_dB)+2 max(Pxx_dB)+2], 'color', 'k', 'LineWidth', 4);
        set(gca, 'ColorOrderIndex', get(gca, 'ColorOrderIndex')-1);
    end
    dpss_validFreqInds= Freq_xx>(stim_freq_shift-freq_half_window_co_Hz) & Freq_xx<(stim_freq_shift+freq_half_window_co_Hz);
    dpss_legalFreqInds= Freq_xx>freqMin;
    outPower= sum(db2pow(Pxx_dB(dpss_validFreqInds)));
    totPower= sum(db2pow(Pxx_dB(dpss_legalFreqInds)));
end