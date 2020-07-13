% function FigS1test_nlIHCfun_PSD(saveFig, LatexDir)
function FigS1test_nlIHCfun_PSD(saveFig, LatexDir)

if ~exist('saveFig', 'var')
    saveFig= 0;
end
if ~exist('LatexDir', 'var')
    LatexDir= ['figures' filesep];
end

figSize_cm= [15 5 13.2 8];
figure_prop_name = {'PaperPositionMode','units','Position', 'Renderer'};
figure_prop_val =  { 'auto'            ,'centimeters', figSize_cm, 'painters'};  % [Xcorner Ycorner Xwidth Ywidth]
figure(1);
clf;
set(gcf,figure_prop_name,figure_prop_val);

plt.xtick_psd_val= [100 500 1e3 2e3 3e3];
plt.xtick_psd_lab= cellfun(@(x) num2str(x), num2cell(plt.xtick_psd_val/1e3), 'uniformoutput', false);

fs= 20e3;
fc= 1e3;
fm= 100;
modDepth= 1;
dur= 1;

[x_sam, t]= helper.create_SAM(fc, fm, fs, modDepth, dur);
t= t*1e3;
x_ihc= 2*max(x_sam)*(normcdf(3*x_sam/max(x_sam))-.5);
lp_nh= helper.get_filter_designfilt('lp', 2e3, fs, 4);
x_ihc= filtfilt(lp_nh, x_ihc);

ax(1)= subplot(221);
plot(t, x_sam, 'Color', helper.get_color('g'));
xlHan_T= xlabel('Time (ms)', 'Units', 'normalized', 'HorizontalAlignment', 'center');
box off;
ttlHan(1)= title('SAM');
ylHan(1)= ylabel('Amplitude', 'Units', 'normalized');

ax(2)= subplot(223);
hold on;
[~,~, lHan(1)]= helper.plot_dpss_psd(x_sam, fs);
[~,~, lHan(2)]= helper.plot_dpss_psd(abs(x_sam), fs);

set(lHan(1), 'color', helper.get_color('g'));
set(lHan(2), 'color', helper.get_color('prp'));
set(lHan, 'linew', 1);

[lgHan, icons]= legend('D(f)', 'S(f)', 'box', 'off', 'Location', 'northwest');
icons(3).XData= mean(icons(3).XData) + [0 +.25];
icons(5).XData= mean(icons(5).XData) + [0 +.25];
lgHan.Position(1:2)= [.15 .33];
set(gca, 'XTick', plt.xtick_psd_val, 'XTickLabel', plt.xtick_psd_lab);

xlHan_F= xlabel('Frequency (kHz)', 'Units', 'normalized', 'HorizontalAlignment', 'center');
ylHan(2)= ylabel('PSD (dB/Hz)', 'Units', 'normalized');

ax(3)= subplot(222);
plot(t, x_ihc, 'color', helper.get_color('g'));
box off;
ttlHan(2)= title('vIHC');

ax(4)= subplot(224);
hold on;
[~,~, lHan(1)]= helper.plot_dpss_psd(x_ihc, fs);
[~,~, lHan(2)]= helper.plot_dpss_psd(abs(x_ihc), fs);

set(lHan(1), 'color', helper.get_color('g'));
set(lHan(2), 'color', helper.get_color('prp'));

set(lHan, 'linew', 1);
xlabel('');
ylabel('');
set(gca, 'XTick', plt.xtick_psd_val, 'XTickLabel', plt.xtick_psd_lab);


linkaxes(ax([1 3]));
xlim(ax(1), [0 2.2/fm]*1e3);
ylim(ax(1), [-2.4 2.4]);
linkaxes(ax([2 4]));
xlim(ax(2), [fm/2 3.5*fc]);
ylim(ax(2), [-60 0]);
set(findall(gcf,'-property','FontSize'),'FontSize', 9);

xlHan_T.Position(1)= 1.1;
xlHan_F.Position(1)= 1.1;

plt.tick_len= [.02 .02];
set(findall(gcf,'-property','TickLength'),'TickLength', plt.tick_len, 'units', 'normalized');
set(ttlHan, 'FontSize', 11);
helper.add_subplot_letter_txp(2, 2, 11, .03, 1.0);

ylHan(1).Position(1)= -.14;
ylHan(2).Position(1)= -.14;

Xcorner_AB= .09;
Xwidth_AB= .41;
Xshift_horz= .08;
Ycorner_AB= .112;
Ywidth_AB= .34;
Yshift_AB= .14;

% B
set(ax(2),'Position',[Xcorner_AB Ycorner_AB Xwidth_AB Ywidth_AB])
drawnow
% A
set(ax(1),'Position',[Xcorner_AB Ycorner_AB+Ywidth_AB+Yshift_AB Xwidth_AB Ywidth_AB])
drawnow
% define new axes for AB
Xwidth_CD=Xwidth_AB;
Ywidth_CD= Ywidth_AB;
Xcorner_CD= Xcorner_AB+Xwidth_AB+Xshift_horz;
Yshift_CD= Yshift_AB;
Ycorner_CD= Ycorner_AB;

% D
set(ax(4),'Position',[Xcorner_CD Ycorner_CD Xwidth_CD Ywidth_CD])
drawnow
% C
set(ax(3),'Position',[Xcorner_CD Ycorner_CD+Ywidth_CD+Yshift_CD Xwidth_CD Ywidth_CD])
drawnow

if saveFig
    saveas(gcf, [LatexDir 'FigS1'], 'epsc');
end