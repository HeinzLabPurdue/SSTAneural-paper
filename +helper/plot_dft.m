% To plot DFT
% function [amp_dB, freq, ax]=plot_dft(vecin, fs, xscale, optional)
% ---------------------------------------------------------------------------
% Example: [Pxx_dB,freq, ax]=plot_dft(vecin, fs, 'log', 'NW', NW...
% 'title', title_str, 'DC', true, 'plot', true, 'nfft', nfft)
% ---------------------------------------------------------------------------
% vecin: input vector, can be either row or column vector.
% fs: Sampling Frequency
% xscale: {'log' (default), 'lin'}
% 'title', title_str
% 'DC', false (default) or true: demean signal if false
% 'phase', false (default) or true: plot phase or not
% 'plot', true (default) or false: do not plot PSD if false
% 'nfft', nfft: for frequency resolution

function varargout=plot_dft(vecin, fs, varargin)

warn_stat= warning('query');
warning('on');

%% parse input
p=inputParser;
x_scale_default= 'log';
y_scale_default= 'log';

WindowType_default= 'rect';
WindowType_valid={'hamming', 'bartlett', 'hann', 'blackman', 'rect', 'rectwin', 'slepian', 'dpss'};
WindowType_check= @(x) any(validatestring(x, WindowType_valid));

title_default= '';
DC_default= false;
plot_default= true;
nfft_default= 2^nextpow2(length(vecin));
phase_default= false;
yRange_default= -100;
sided_default= 1;
xunit_default= 'hz';

addRequired(p, 'vecin', @isnumeric);
addRequired(p, 'fs', @isnumeric);
addOptional(p,'window', WindowType_default, WindowType_check);
addParameter(p,'xscale', x_scale_default, @ischar);
addParameter(p,'yscale', y_scale_default, @ischar);
addParameter(p,'title', title_default, @ischar);
addParameter(p,'DC', DC_default, @islogical);
addParameter(p,'phase', phase_default, @islogical);
addParameter(p,'plot', plot_default, @islogical);
addParameter(p,'NFFT', nfft_default, @isnumeric);
addParameter(p,'yRange', yRange_default, @isnumeric);
addParameter(p,'sided', sided_default, @isnumeric);
addParameter(p,'xunit', xunit_default, @ischar);

p.KeepUnmatched = true;
parse(p,vecin,fs,varargin{:})

%%
vecin= vecin(:);
L = length(vecin);

if ~p.Results.DC
    dc_val= mean(vecin);
    vecin=vecin-dc_val;
end

switch lower(p.Results.window)
    case {'slepian', 'dpss'}
        t_window=dpss(L, 1, 1);
    case 'hamming'
        t_window=hamming(L);
    case 'bartlett'
        t_window=bartlett(L);
    case 'hann'
        t_window=hann(L);
    case 'blackman'
        t_window=blackman(L);
    case {'rect', 'rectwin'}
        t_window=rectwin(L);
end


if strcmpi(p.Results.yscale, 'dbspl')
    if ~any(ismember(lower(p.UsingDefaults), {'nfft'}))
        warning('Not using nfft provided by user as yscale is dbspl. Instead using nfft=L.');
    end
    nfft= L;
    if ~ismember(p.Results.window, {'rect', 'rectwin'})
        warning('Not using rectangular window for dbspl. dBSPL values may not be exact. ');
    end
else
    nfft= p.Results.NFFT;
end

E_window= rms(t_window);
Y= fft(vecin.*t_window, nfft)/E_window;
P2= Y/L;

if p.Results.sided==1
    amp= abs(P2(1:ceil(nfft/2+1)));
    amp(2:end-1)= 2*amp(2:end-1);
    amp_dB= db(amp); % 20*log10(amp) = Power
    freq= linspace(0,fs/2,length(amp_dB));
    xscale= p.Results.xscale;
    xl_val= [fs/2/nfft fs/2];
elseif p.Results.sided==2
    amp= fftshift(abs(P2));
    amp_dB= db(amp); % 20*log10(amp) = Power
    if rem(length(amp_dB), 2)==1
        freq =linspace(-fs/2,fs/2,length(amp_dB));
    else 
        freq= linspace(0,fs/2,length(amp_dB)/2+1);
        freq= [-fliplr(freq(2:end-1)), freq];
    end
    xscale= 'lin';
    % Warn if not using default window (means user wanted log), but force
    % linear since two-sided
    if ~strcmp(xscale, p.Results.xscale) && ~any(ismember(p.UsingDefaults, 'xscale'))
        warning('Using linear x-scale insteald of %s as two-sided FFT', p.Results.xscale);
    end
    xl_val= [-fs/2 fs/2];
end

if ismember(p.Results.yscale, {'log', 'db'})
    amp_to_plot= amp_dB;
    ylab_str= '|P1(f)| dB';
elseif ismember(p.Results.yscale, {'lin', 'mag'})
    amp_to_plot= amp;
    ylab_str= '|P1(f)| mag';
elseif ismember(p.Results.yscale, {'dbspl'})
    amp_to_plot= helper.dbspl(amp/sqrt(2)); % RMS = AMP/sqrt(2);
    ylab_str= '|P1(f)| in dB SPL';
end

if p.Results.plot
    if p.Results.phase
        ax(1)=subplot(211);
    end
    if ismember(p.Results.xunit, {'khz', 'k'})
        freq= freq/1e3;
        xlim_div= 1e3;
        xlab_str= 'Frequency (kHz)';
    elseif ismember(p.Results.xunit, {'Hz', 'hz'})
        xlim_div= 1;
        xlab_str= 'Frequency (Hz)';
    end
    
    lHan=plot(freq, amp_to_plot, 'linew', 2);
%     grid on;
    xlabel(xlab_str);
    ylabel(ylab_str);
    title(p.Results.title);
    
    if p.Results.phase
        ax(2)=subplot(212);
        semilogx(freq, unwrap(angle(P2(1:ceil(nfft/2+1)))), 'linew', 2);
        title('Phase Plot');
        linkaxes(ax, 'x');
        grid on;
    end
    
    set(gca, 'xscale', xscale);
    xlim(xl_val/xlim_div);
    if ismember(p.Results.yscale, {'log', 'db', 'dbspl'}) &&  p.Results.yRange>0
        ylim([max(amp_to_plot)-p.Results.yRange+10 max(amp_to_plot)+10]);
    end
else 
    lHan= nan;
end

if nargout
    varargout{1}= amp_to_plot;
    varargout{2}= freq;
    varargout{3}= lHan;
end

warning(warn_stat(1).state);
