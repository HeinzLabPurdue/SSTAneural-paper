% function FigS1_graphical_illustration(saveFig, LatexDir)
function FigS1_graphical_illustration(saveFig, LatexDir)

if ~exist('saveFig', 'var')
    saveFig= 0;
end
if ~exist('LatexDir', 'var')
    LatexDir= ['figures' filesep];
end

fc= 100;
fm= 20;
fs= 1e3;
[x, t]= create_SAM(fc, fm, fs, 1, 1);
t_ms= t*1e3;

p_t= x .* (x>0);
n_t= -x .* (-x>0);
s_t= (p_t+n_t)/2;
d_t= (p_t-n_t)/2;
a_t= rms(d_t)*hilbert(d_t);
e_t= abs(a_t)/sqrt(2);
phi_t= rms(d_t) * cos(angle(a_t));


xtick_vals_freq= [20 100 200 400];
xtick_labs_freq= cellfun(@(x) num2str(x), num2cell(xtick_vals_freq), 'UniformOutput', false);
lw= 1.5;
fSize= 12;

figSize_cm= [15 5 19.05 13];
figure_prop_name = {'PaperPositionMode', 'units', 'Position', 'Renderer'};
figure_prop_val =  { 'auto', 'centimeters', figSize_cm, 'painters'};  % [Xcorner Ycorner Xwidth Ywidth]
figure(1);
clf;
set(gcf,figure_prop_name,figure_prop_val);

sp_ax(1)= subplot(321);
hold on;
plot(t_ms, p_t, 'linew', lw);
plot(t_ms, n_t, '-.', 'linew', lw);
title('Time waveform');

[lgHan, icons]= legend('p(t)', 'n(t)', 'box', 'off', 'FontSize', fSize);
icons(3).XData= mean(icons(3).XData) + [0 +.25];
icons(5).XData= mean(icons(5).XData) + [0 +.25];
lgHan.Position(1:2)= [.37 .85];


sp_bx(1)= subplot(322);
hold on;
[~,~,lHan1]= plot_dpss_psd(p_t, fs);
[~,~,lHan2]= plot_dpss_psd(n_t, fs);
set(lHan1, 'linestyle', '-', 'linew', lw);
set(lHan2, 'linestyle', '-.', 'linew', lw);
title('Spectrum');
ylabel('');
xlabel('');
set(gca, 'xtick', xtick_vals_freq, 'XTickLabel', xtick_labs_freq);

[lgHan, icons]= legend('P(f)', 'N(f)', 'box', 'off', 'FontSize', fSize);
icons(3).XData= mean(icons(3).XData) + [0 +.25];
icons(5).XData= mean(icons(5).XData) + [0 +.25];
lgHan.Position(1:2)= [.87 .85];

sp_ax(2)= subplot(323);
hold on;
plot(t_ms, s_t, 'linew', lw);
plot(t_ms, d_t, '-.', 'linew', lw);
ylabel('Amplitude');

[lgHan, icons]= legend('s(t)', 'd(t)', 'box', 'off', 'FontSize', fSize);
icons(3).XData= mean(icons(3).XData) + [0 +.25];
icons(5).XData= mean(icons(5).XData) + [0 +.25];
lgHan.Position(1:2)= [.37 .54];

sp_bx(2)= subplot(324);
hold on;
[~,~,lHan1]= plot_dpss_psd(s_t, fs);
[~,~,lHan2]= plot_dpss_psd(d_t, fs);
set(lHan1, 'linestyle', '-', 'linew', lw);
set(lHan2, 'linestyle', '-.', 'linew', lw);
xlabel('');
set(gca, 'xtick', xtick_vals_freq, 'XTickLabel', xtick_labs_freq);

[lgHan, icons]= legend('S(f)', 'D(f)', 'box', 'off', 'FontSize', fSize);
icons(3).XData= mean(icons(3).XData) + [0 +.25];
icons(5).XData= mean(icons(5).XData) + [0 +.25];
lgHan.Position(1:2)= [.87 .54];

sp_ax(3)= subplot(325);
hold on;
plot(t_ms, e_t, 'linew', lw);
plot(t_ms, phi_t, '-.', 'linew', lw);
xlabel('Time (ms)');

[lgHan, icons]= legend('e(t)', '\phi(t)', 'box', 'off', 'FontSize', fSize);
icons(3).XData= mean(icons(3).XData) + [0 +.25];
icons(5).XData= mean(icons(5).XData) + [0 +.25];
lgHan.Position(1:2)= [.37 .24];

sp_bx(3)= subplot(326);
hold on;
[~,~,lHan1]= plot_dpss_psd(e_t, fs);
[~,~,lHan2]= plot_dpss_psd(phi_t, fs);
set(lHan1, 'linestyle', '-', 'linew', lw);
set(lHan2, 'linestyle', '-.', 'linew', lw);
ylabel('');
ylim([-50 0]);
set(gca, 'xtick', xtick_vals_freq, 'XTickLabel', xtick_labs_freq);

[lgHan, icons]= legend('E(f)', '\Phi(f)', 'box', 'off', 'FontSize', fSize);
icons(3).XData= mean(icons(3).XData) + [0 +.25];
icons(5).XData= mean(icons(5).XData) + [0 +.25];
lgHan.Position(1:2)= [.87 .24];

linkaxes(sp_ax);
xlim(sp_ax(1), [0 2.2/fm*1e3]);
set(findall(gcf,'-property','FontSize'),'FontSize', fSize);


linkaxes(sp_bx);
xlim(sp_bx(1), [10 fs/2]);
ylim(sp_bx(1), [-60 -10]);

add_subplot_letter(3, 2, 14);

%% define new axes for AB
Xcorner_X= .075;
Xwidth_X= .39;
Xshift_X= .12;
Ycorner_X= .095;
Ywidth_X= .24;
Yshift_X= .07;

% E
set(sp_ax(3),'Position',[Xcorner_X, Ycorner_X, Xwidth_X Ywidth_X])
drawnow
% C
set(sp_ax(2),'Position',[Xcorner_X, Ycorner_X+Ywidth_X+Yshift_X, Xwidth_X, Ywidth_X])
drawnow
% A
set(sp_ax(1),'Position',[Xcorner_X, Ycorner_X+2*Ywidth_X+2*Yshift_X, Xwidth_X, Ywidth_X])
drawnow
% F
set(sp_bx(3),'Position',[Xcorner_X+Xwidth_X+Xshift_X, Ycorner_X, Xwidth_X, Ywidth_X])
drawnow
% D
set(sp_bx(2),'Position',[Xcorner_X+Xwidth_X+Xshift_X, Ycorner_X+Ywidth_X+Yshift_X, Xwidth_X, Ywidth_X])
drawnow
% B
set(sp_bx(1),'Position',[Xcorner_X+Xwidth_X+Xshift_X, Ycorner_X+2*Ywidth_X+2*Yshift_X, Xwidth_X, Ywidth_X])
drawnow

if saveFig 
    saveas(gcf, [LatexDir 'FigS1'], 'epsc');
end