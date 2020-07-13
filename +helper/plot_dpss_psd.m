% To plot PSD estimated using DPSS sequences.
% function [Pxx_dB,freq, ax]=plot_dpss_psd(data_in, fs, optional)
% ---------------------------------------------------------------------------
% Example: [Pxx_dB,freq, ax]=plot_welch_psd(data_in, fs, 'NW', NW, 'nfft', nfft,  ...
% 'plot', true, 'plotconf', true, 'xscale', 'log', 'yscale', 'dB', 'title', title_str...
% 'xunit', 'khz', 'DC', true, 'yrange', 40)
% ---------------------------------------------------------------------------
% data_in: input vector, can be either row or column vector.
% fs: Sampling Frequency
% 'xscale', {'log' || 'lin'}
% 'yscale', {'dB', 'log' || 'mag', 'lin'}
% 'title', title_str
% 'DC', {true or false}: demean signal if false (default false)
% 'plot', {true or false}: do not plot PSD if false (default true)
% 'plotconf', {true or false}: do not plot PSD Confidence Interval if false (default false)
% 'nfft', nfft: for frequency resolution
% 'NW', NW: for time-bandwidth half-product
% 'xunit', {'Hz', 'hz' || 'khz', 'k'}: kHz or Hz for frequency axis
% 'yrange', 40: y-range for PSD (default: no limit)

function varargout=plot_dpss_psd(data_in, fs, varargin)

%% parse input
p=inputParser;
NW_default= 2.5;

xscale_default= 'log';
yscale_default= 'dB';
title_default= '';
DC_default= false;
plotConf_default=false;
plot_default= true;
nfft_default= 2^nextpow2(length(data_in));
yrange_default= -1;
ConfidenceLevel_default=.95;
xunit_default= 'hz';
norm_default= false;

addRequired(p, 'data_in', @isnumeric);
addRequired(p, 'fs', @isnumeric);
addParameter(p,'NW', NW_default, @isnumeric)
addParameter(p,'NFFT', nfft_default, @isnumeric);
addParameter(p,'xscale', xscale_default, @ischar);
addParameter(p,'yscale', yscale_default, @ischar);
addParameter(p,'title', title_default, @ischar);
addParameter(p,'DC', DC_default, @islogical);
addParameter(p,'norm', norm_default, @islogical);
addParameter(p,'plot', plot_default, @islogical);
addParameter(p,'plotConf', plotConf_default, @islogical);
addParameter(p,'yrange', yrange_default, @isnumeric);
addParameter(p,'ConfidenceLevel', ConfidenceLevel_default, @isnumeric);
addParameter(p,'xunit', xunit_default, @ischar);

p.KeepUnmatched = true;
parse(p,data_in,fs,varargin{:})


%% actual code
if isrow(data_in)
    data_in= data_in';
end

if ~p.Results.DC
    data_in=data_in-mean(data_in,1);
end

[Pxx_pow, freq, Pxx_ci_lin] = pmtm(data_in, p.Results.NW, p.Results.NFFT, fs, 'ConfidenceLevel', p.Results.ConfidenceLevel);

if p.Results.norm
    % Normalize to get PSD
    Pxx_pow= Pxx_pow/p.Results.NFFT*fs;
    Pxx_ci_lin= Pxx_ci_lin/p.Results.NFFT*fs;
    % else
    %     Pxx_lin= Pxx_lin; % /length(data_in). Shouldn't change pmtm output.
end

freq_inds2use= freq~=0;
freq=freq(freq_inds2use);
Pxx_pow=Pxx_pow(freq_inds2use, :);
Pxx_ci_lin=Pxx_ci_lin(freq_inds2use,:);

if ismember(p.Results.xunit, {'khz', 'k'})
    freq= freq/1e3;
    xlim_div= 1e3;
    xlab_str= 'Frequency (kHz)';
elseif ismember(p.Results.xunit, {'Hz', 'hz'})
    xlim_div= 1;
    xlab_str= 'Frequency (Hz)';
end

Pxx_dB= pow2db(Pxx_pow); % Previously: Pxx_lin is already A^2. So 10log10(*) for power => db(*)/2.
Pxx_ci_dB= pow2db(Pxx_ci_lin);
yl_val = nan;

%% Plot here
if p.Results.plot
    if ismember(lower(p.Results.yscale), {'log', 'db'})
        Pxx_to_plot= Pxx_dB;
        Pxx_ci_to_plot= Pxx_ci_dB;
        ylab_str = 'PSD (dB/Hz)';
        
        if p.Results.yrange>0
            yl_val = [max(Pxx_to_plot)-p.Results.yrange max(Pxx_to_plot)+5];
        end
        
    elseif ismember(lower(p.Results.yscale), {'lin', 'mag'})
        Pxx_to_plot= Pxx_pow;
        Pxx_ci_to_plot= Pxx_ci_lin;
        ylab_str = 'PSD (Amp_l_i_n/Hz)';
        
        if p.Results.yrange>0 % means not default value
            warning('The option ''yrange'' is valid only for log y-scale');
        end
    elseif ismember(lower(p.Results.yscale), {'dbspl', 'spl'})
        Pxx_to_plot= dbspl(sqrt(Pxx_pow)/sqrt(2));
        Pxx_ci_to_plot= dbspl(sqrt(Pxx_ci_lin)/sqrt(2));
        ylab_str = 'Power in dB SPL';
        
        if p.Results.yrange>0
            yl_val = [max(Pxx_to_plot)-p.Results.yrange max(Pxx_to_plot)+5];
        end
    end
    
    % Plot PSD and CI
    ax=plot(freq, Pxx_to_plot, '-', 'linew', 2);
    
    if nargout==4 || p.Results.plotConf
        hold on;
        ax= [ax; ax];
        for plotVar=1:size(Pxx_to_plot,2)
            bx= fill([freq;flipud(freq)], [Pxx_ci_to_plot(:,2*(plotVar-1)+1); flipud(Pxx_ci_to_plot(:,2*plotVar))], get(ax(plotVar), 'color'));
            bx.EdgeAlpha= 0;
            bx.FaceAlpha=.4;
            set(gca, 'ColorOrderIndex', max(1, get(gca, 'ColorOrderIndex')-1));
            ax(size(Pxx_to_plot,2) + plotVar)= bx;
        end
        %         hold off;
    end
    
    % Plotting parameters
    set(gca, 'xscale', p.Results.xscale);
    xlabel(xlab_str);
    ylabel(ylab_str);
    title(p.Results.title);
    xlim([fs/p.Results.NFFT fs/2]/xlim_div);
    %     grid on;
    if ~isnan(yl_val)
        ylim(yl_val);
    end
    
else
    ax= nan;
end

%% Assign outputs
if nargout
    varargout{1}= Pxx_dB;
    varargout{2}= freq;
    varargout{3}= ax;
    varargout{4}= Pxx_ci_dB;
end