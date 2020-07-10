function PSDout= plot_snap_fft_natural_vowel...
    (uRatePos, uRateNeg, fsResp, stim, fsStim, tStart, tEnd, tc_data, doPlot)

figSize_cm= [15 5 13.2 8];
figure_prop_name = {'PaperPositionMode', 'units', 'Position', 'Renderer'};
figure_prop_val =  { 'auto', 'centimeters', figSize_cm, 'painters'};  % [Xcorner Ycorner Xwidth Ywidth]
figure(1);
clf;
set(gcf,figure_prop_name, figure_prop_val);


tcLIMs= [10 75];
nStart=max(1, round(tStart*fsResp));
nEnd=min(length(uRatePos), round(tEnd*fsResp));
uRate_Pos_Snap=uRatePos(nStart:nEnd);
uRate_Neg_Snap=uRateNeg(nStart:nEnd);
t_resp=(1:length(uRatePos))/fsResp;
tSnap=(nStart:nEnd)/fsResp;

tc_data.TCfit= helper.trifilt(tc_data.TCfit, 5);

xLimVals=[75 10e3]/1e3;
xTicks= [.1 1 max(xLimVals)];
xTickLabel= cellfun(@(x) num2str(x), num2cell(xTicks), 'uniformoutput', false);

uRateTFS=(uRatePos-uRateNeg)/2;
uRateTFS_snap=uRateTFS(nStart:nEnd);

uRateSUM=(uRatePos+uRateNeg)/2;
uRateSUM_snap=uRateSUM(nStart:nEnd);

t_stim_snap=(max(1, round(fsStim*tStart)):round(fsStim*tEnd))/fsStim;
stim_snap=stim(max(1, round(fsStim*tStart)):round(fsStim*tEnd));

lw=.75;
lw3= 1.2;

nfft2=2^nextpow2(length(uRatePos)); % if 2056, decent and can save epsc file
nfft3= 2^nextpow2(length(uRatePos));
NW= 1.5;
keep_mean_in_psd=false;
x_scale_PSD= 'log';
y_scale_PSD= 'log';

if tStart==.48 && tEnd==.58
    F0_freq= 110/1e3; % true value = 130, adjust for plotting
    F1_freq= 220/1e3; % true value = 190, adjust for plotting
    F2_freq= 1320/1e3;
    textYval= -20;
else
    error('Formant values are hardcoded for 0.48 - 0.58 s');
end

%%
if doPlot
    co=get(gca, 'colororder');
    ax(1)= subplot(2,3,1:2);
    tStim=(1:length(stim))/fsStim;
    Astim=.7*max([uRatePos, uRateNeg]);
    
    hold on;
    plot(tStim, Astim*stim-Astim, '-', 'color', 2*helper.get_color('gray'), 'linew', lw);
    lg(3)= plot(t_stim_snap, Astim*stim_snap-Astim, '-', 'color', helper.get_color('k'), 'linew', lw);
    plot(t_resp, uRatePos, '-', 'color', 2*helper.get_color('gray'), 'linew', lw);
    plot(t_resp, -uRateNeg, '-', 'color', 2*helper.get_color('gray'), 'linew', lw);
    lg(1)= plot(tSnap, uRate_Pos_Snap, '-', 'color', helper.get_color('b'), 'linew', lw);
    lg(2)= plot(tSnap, -uRate_Neg_Snap, '-', 'color', helper.get_color('r'), 'linew', lw);
    
    [lgHan, icons]= legend(lg, 'p(t)', '-n(t)', 'Stim', 'location', 'southeast', 'box', 'off', 'interpreter', 'tex');
    lgHan.Box= 'off';
    lgHan.Position(1:2)= [.45 .6];
    
    icons(4).XData= mean(icons(4).XData) + [0.1 .3];
    icons(4).LineWidth= 1;
    icons(6).XData= mean(icons(6).XData) + [0.1 .3];
    icons(6).LineWidth= 1;
    icons(8).XData= mean(icons(8).XData) + [0.1 .3];
    icons(8).LineWidth= 1;

    
    ylabel('Amplitude');
    box off;
    xlabel('Time (s)');
    xlim([max(.7*tStart, tStart-.02) .01+min([1.4*tEnd, tEnd+.04, tStim(end)+.03])]);
    ylim([-1.2*Astim .5*Astim])
    txtHan(1)= text(.05, 1.075, 'A. Stimulus & Response', 'interpreter', 'tex', 'units', 'normalized');
    
    %%
    ax(2)= subplot(2,3,3);
    yyaxis left;
    yRange_pdf= 50;
    [Pxx_stim, ~, lHan]= helper.plot_dft(helper.gen_rescale(stim_snap, 70), fsStim, 'xscale', x_scale_PSD, 'DC', keep_mean_in_psd, 'xunit', 'khz', 'yrange', yRange_pdf, 'yscale', 'dbspl');
    hold on;
    text(F0_freq, max(Pxx_stim)+5, 'F_0', 'HorizontalAlignment', 'center');
    text(F1_freq, max(Pxx_stim)+2, 'F_1', 'HorizontalAlignment', 'center');
    text(F2_freq, max(Pxx_stim)+3-15, 'F_2', 'HorizontalAlignment', 'center');
    ylim([max(Pxx_stim)-yRange_pdf max(Pxx_stim)+10]);
    set(lHan, 'linew', lw3);
    
    set(gca, 'xtick', xTicks, 'xticklabel', xTickLabel);
    title('');
    txtHan(2)= text(.05, 1.075, 'B. Stimulus', 'interpreter', 'tex', 'units', 'normalized');
    xlim(xLimVals)
    grid off;
    ylabel('DFT (dB SPL)');
    xlabel('');
    set(gca, 'xtick', xTicks, 'xticklabel', xTickLabel);
    set(gca, 'ycolor', 'k');
    
    yyaxis right;
    plot(tc_data.freqkHz, tc_data.TCfit, 'k', 'linew', lw3)
    set(gca, 'ycolor', 'k')
    ylim(tcLIMs);
    box off;
    
    %%
    ax(3)=subplot(2,3,4);
    
    [Pxx_lin_pos,~, lPlot2]=helper.plot_dpss_psd(uRate_Pos_Snap, fsResp, 'xscale', x_scale_PSD, 'NW', NW, 'nfft', nfft2, 'DC', keep_mean_in_psd, 'yscale', y_scale_PSD, 'xunit', 'khz');
    hold on;
    text(F0_freq, textYval-5, 'F_0', 'HorizontalAlignment', 'center');
    text(F2_freq, textYval-5, 'F_2', 'HorizontalAlignment', 'center');
    
    title('');
    txtHan(3)= text(.05, 1.05, 'C. P(f)', 'interpreter', 'tex', 'units', 'normalized');
    xlim(xLimVals)
    set(gca, 'xtick', xTicks, 'xticklabel', xTickLabel);
    grid off;
    
    yyaxis right;
    plot(tc_data.freqkHz, tc_data.TCfit, 'k', 'linew', lw3)
    set(gca, 'ycolor', 'k', 'yticklabel', '')
    xlabel('');
    ylim(tcLIMs);
    
    yyaxis left;
    set(lPlot2, 'color', co(1,:), 'linew', lw3);
    set(gca, 'ycolor', 'k');
    box off;
    
    ax(4)=subplot(2,3,5);
end

hold on;
[PSDout.Pxx_tfs,PSDout.freq_tfs, lPlot_comp]=helper.plot_dpss_psd(uRateTFS_snap, fsResp, 'xscale', x_scale_PSD, 'NW', NW, ...
    'nfft', nfft2, 'DC', keep_mean_in_psd, 'yscale', y_scale_PSD, 'plot', doPlot==true, 'xunit', 'khz');
[PSDout.Pxx_sum,~, lPlot_sum]=helper.plot_dpss_psd(uRateSUM_snap, fsResp, 'xscale', x_scale_PSD, 'NW', NW, ...
    'nfft', nfft2, 'DC', keep_mean_in_psd, 'yscale', y_scale_PSD, 'plot', doPlot==true, 'xunit', 'khz');
PSDout.nSpikes= (sum(uRate_Pos_Snap)+sum(uRate_Neg_Snap))/2;
set(lPlot_comp, 'color', helper.get_color('g'), 'linew', lw3);
set(lPlot_sum, 'color', helper.get_color('prp'), 'linew', lw3);

ylabel('');

if doPlot
    hold on;
    text(F0_freq, textYval-5, 'F_0', 'HorizontalAlignment', 'center');
    text(F2_freq, textYval-5, 'F_2', 'HorizontalAlignment', 'center');
    set(gca, 'xtick', xTicks, 'xticklabel', xTickLabel, 'YTickLabel', '');
    xlabel('Frequency (kHz)');
    %     title('');
    txtHan(4)= text(.05, 1.05, 'D. D(f) & S(f)', 'interpreter', 'tex', 'units', 'normalized');
    xlim(xLimVals)
    grid off;
    
    yyaxis right;
    plot(tc_data.freqkHz, tc_data.TCfit, 'k', 'linew', lw3)
    set(gca, 'ycolor', 'k', 'yticklabel', '')
    ylim(tcLIMs);
    
    yyaxis left;
    box off;
    [legTxtHan, icons]= legend([lPlot_comp lPlot_sum], 'D(f)', 'S(f)', 'Location', 'northeast', 'box', 'off');
    legTxtHan(1).Position(1)= .48;
    icons(3).XData= mean(icons(3).XData) + [0.1 0.3];
    icons(5).XData= mean(icons(5).XData) + [0.1 0.3];

    %%
    ax(5)= subplot(2,3,6);
    
    hold on;
    [~, ~, ll]=helper.plot_dpss_psd(uRateTFS_snap./abs(hilbert(uRateTFS_snap)), fsResp, 'xscale', x_scale_PSD, 'NW', NW, 'nfft', nfft3, 'DC', keep_mean_in_psd, 'yscale', y_scale_PSD, 'plot', true, 'xunit', 'khz');
    set(ll', 'color', helper.get_color('g'), 'DisplayName', '\Phi(f)', 'linew', lw3)
    hold on;
    [~, ~, lPlot4]=helper.plot_dpss_psd(abs(hilbert(uRateTFS_snap)), fsResp, 'xscale', x_scale_PSD, 'NW', NW, 'nfft', nfft3, 'DC', keep_mean_in_psd, 'yscale', y_scale_PSD, 'plot', true, 'xunit', 'khz');
    set(lPlot4, 'color', helper.get_color('prp'), 'DisplayName', 'E(f)', 'linew', lw3);
    text(F0_freq, textYval-4, 'F_0', 'HorizontalAlignment', 'center');
    text(F2_freq, textYval-5, 'F_2', 'HorizontalAlignment', 'center');
    
    set(gca, 'xtick', xTicks, 'xticklabel', xTickLabel, 'YTickLabel', '');
    xlabel('');
    ylabel('');
    title('');
    txtHan(5)= text(.05, 1.05, 'E. \Phi(f) & E(f)', 'interpreter', 'tex', 'units', 'normalized');
    grid off;
    box off;
    
    yyaxis right;
    plot(tc_data.freqkHz, tc_data.TCfit, '-k', 'linew', lw3, 'DisplayName', 'TC')
    set(gca, 'ycolor', 'k')
    ylim(tcLIMs);
    xlim(xLimVals)
    ylabel('TC \theta (dB SPL)', 'interpreter', 'tex');
    
    [legTxtHan, icons]= legend('show', 'location', 'northeast', 'box', 'off', 'interpreter', 'tex');
    legTxtHan(1).Position(1)= .77;
    icons(4).XData= mean(icons(4).XData) + [0.1 0.3];
    icons(6).XData= mean(icons(6).XData) + [0.1 0.3];
    icons(8).XData= mean(icons(8).XData) + [0.1 0.3];

    yyaxis left;
    
    linkaxes(ax(2:4), 'x');
    linkaxes(ax(3:5), 'y');
    
    yyaxis left;
    ylim(ax(5), [max(Pxx_lin_pos)-30 max(Pxx_lin_pos)+10]);
    
    
    %% define new axes for AB
    Xshift_horz= .06;
    Xwidth_C= .235;
    Ywidth_C= .34;
    Xcorner_C= .09;
    Yshift_C= .15;
    Ycorner_C= .113;
    
    % A
    set(ax(1),'Position',[Xcorner_C, Ycorner_C+Ywidth_C+Yshift_C, 2*Xwidth_C+Xshift_horz-.02, Ywidth_C])
    drawnow

    % B
    set(ax(2),'Position',[Xcorner_C+2*Xwidth_C+2*Xshift_horz, Ycorner_C+Ywidth_C+Yshift_C, Xwidth_C, Ywidth_C])
    drawnow

    % C
    set(ax(3),'Position',[Xcorner_C, Ycorner_C, Xwidth_C, Ywidth_C])
    drawnow

    % D
    set(ax(4),'Position',[Xcorner_C+Xwidth_C+Xshift_horz, Ycorner_C, Xwidth_C, Ywidth_C])
    drawnow

    % E
    set(ax(5),'Position',[Xcorner_C+2*Xwidth_C+2*Xshift_horz, Ycorner_C, Xwidth_C, Ywidth_C])
    drawnow
    
    set(findall(gcf,'-property','FontSize'),'FontSize', 9);
    set(txtHan,'FontSize', 11);
    plt.tick_len= [.025 .025];
    set(findall(gcf,'-property','TickLength'),'TickLength', plt.tick_len, 'units', 'normalized');
end