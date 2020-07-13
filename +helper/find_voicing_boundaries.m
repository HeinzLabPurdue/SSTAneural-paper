function voiced_boundaries = find_voicing_boundaries(sig, fs, doPlot, ENV_thresh)

if ~exist('doPlot', 'var')
    doPlot = false;
end
if ~exist('ENV_thresh', 'var')
    ENV_thresh= .1;
end
fSize = 16;

sig= sig/rms(sig);
tSig=((1:length(sig))/fs)';

f_lp_co= 250;
[b,a]= butter(2, f_lp_co/(fs/2), 'low');
[b_env, a_env]= butter(2, 32/(fs/2), 'low');

y_lp = filtfilt(b,a,sig);
y_lp_env= filtfilt(b_env, a_env, y_lp.*(y_lp>0));

% [nCounts, nCenter]= histcounts(y_lp_env, 200);
% ENV_thresh= nCenter(find(cumsum(nCounts./sum(nCounts)) >=.3, 1));

voicing_mask= ones(size(tSig)) .* (y_lp_env>ENV_thresh);

rising_edges  = find(voicing_mask(2:end)>voicing_mask(1:end-1));
falling_edges = find(voicing_mask(2:end)<voicing_mask(1:end-1));
if numel(falling_edges)>numel(rising_edges)
    falling_edges= falling_edges(2:end);
end

voiced_boundaries = [rising_edges(:), falling_edges(:)]/fs;

if doPlot
    figure(1);
    clf;
    
    subplot(211)
    plot(tSig, sig);
    ylabel('Signal');
    
    yyaxis right;
    l= plot(tSig, voicing_mask, 'LineWidth', 3);
    ylim([-.5 1.5]);
    ylabel('Voicing mask');
    set(gca, 'YColor', get(l, 'color'), 'FontSize', fSize);
    
    subplot(212)
    plot(tSig, y_lp);
    ylabel(sprintf('Lowpassed at %.0f Hz', f_lp_co));
    
    yyaxis right;
    plot(tSig, y_lp_env, 'LineWidth', 3);
    hold on;
    plot(tSig, ENV_thresh*ones(size(tSig)), 'k--', 'LineWidth', 3);
    ylabel('HWR+LP + thresholding');
    
    
    set(gca, 'FontSize', fSize);
    
    xlabel('time (sec)')
    
    set(gcf, 'Units', 'inches', 'Position', [1 1 10 7]);
%     saveas(gcf, 'voicing_detect', 'png');
end