% function Fig11_demodulation_demo(saveFig, LatexDir)
function Fig11_demodulation_demo(saveFig, LatexDir)

if ~exist('saveFig', 'var')
    saveFig= 0;
end
if ~exist('LatexDir', 'var')
    LatexDir= ['figures' filesep];
end


%% Create a FM
% define stim params
stim.fs= 10e3;
stim.dur= 2; %188e-3;
stim.siglen= round(stim.fs*stim.dur);
stim.freqStart= 400;
stim.freqEnd= 800;
stim.t= (1:stim.siglen)/stim.fs;
stim.phi_t= stim.freqStart.*stim.t + (stim.freqEnd-stim.freqStart)/2*(stim.t).^2/stim.dur; % Formula for phase in terms of frequency
sig_org= sin(2 * pi * stim.phi_t) + helper.create_sinusoid([1400 2000], stim.fs, stim.dur)';
stim.freqTrajectory= linspace(stim.freqStart, stim.freqEnd, stim.siglen);

anl.tWindow= 100e-3;
anl.nw= 1.5;
anl.dft_sided= 2;
anl.yrange= 100;
anl.dc= true;

%% Demodulate ORG signal numerically (i.e. Phi = int_ (Freq . dt)) where Freq can be any trajectories
mod_empirical= -cumtrapz(stim.freqTrajectory)/stim.fs;
sig_demod_empirical= sig_org .* exp(2*pi*1j*mod_empirical); % sqrt_2 times to conserve power after hilbert transform

y= sig_demod_empirical; %.* exp(2*pi*1j*stim.freqStart.*stim.t);

[A_yy, Freq_yy, ll]= helper.plot_dft(y, stim.fs, 'sided', anl.dft_sided, 'dc', anl.dc, 'yrange', anl.yrange, 'plot', false);

%% Test function
freq_window_Hz= 1;
use_DFT0_DPSS1= 0;
[outPower_fun, totPower_fun]= helper.get_freq_trajectory_power(sig_org, stim.fs, stim.freqTrajectory, false, freq_window_Hz, 100, [], use_DFT0_DPSS1);

validFreqInds= Freq_yy>-freq_window_Hz & Freq_yy<freq_window_Hz;
outPower_psd= sum(db2mag(A_yy(validFreqInds)).^2*numel(y)/numel(A_yy));
totPower_psd= sum(db2mag(A_yy).^2*numel(y)/numel(A_yy));

fprintf('-----> FUN: traj=%.1f, TOT=%.1f \n', outPower_fun, totPower_fun);
fprintf('-----> PSD: traj=%.1f, TOT=%.1f \n', outPower_psd, totPower_psd);

%%
anl.tWindow= 100e-3;
anl.fracOVlap= .99;
anl.nfft= 2^(3+nextpow2(anl.tWindow*stim.fs));

lw= .8;
lw2= 1.5;
yl_val= [-3.0 3.0];
ytick_vals= -2:1:2; % already in kHz
ytick_labs= cellfun(@(x) num2str(x), num2cell(ytick_vals), 'UniformOutput', false);


figSize_cm= [15 5 13.2 8];
figure_prop_name = {'PaperPositionMode','units','Position', 'Renderer'};
figure_prop_val =  { 'auto'            ,'centimeters', figSize_cm, 'painters'};  % [Xcorner Ycorner Xwidth Ywidth]
figure(11);
clf;
set(gcf,figure_prop_name,figure_prop_val);


ax(1)= subplot(221);
spectrogram(sig_org+1j*eps, blackman(anl.tWindow *stim.fs), round(anl.tWindow*anl.tWindow*anl.fracOVlap), anl.nfft, 'centered', 'yaxis', stim.fs);
hold on;
lformant_Han(1)= line(stim.t, stim.freqTrajectory/1e3+.2, 'LineWidth', lw);
lformant_Han(2)= line(stim.t, stim.freqTrajectory/1e3-.2, 'LineWidth', lw);
set(lformant_Han, 'color', helper.get_color('r'), 'LineStyle', '--');
ylim(yl_val);
colorbar off;
xlabel('');
set(gca, 'YTick', ytick_vals, 'YTickLabel', ytick_labs);
ylabel('')

ax(2)= subplot(222);
[Px, Fx]= helper.plot_dft(sig_org, stim.fs, 'sided', 2, 'dc', true, 'plot', false);
xlabel('');
plot(Px, Fx/1e3, 'LineWidth', lw);
hold on;
plot([-30 -30], [min(stim.freqTrajectory)/1e3-.15 max(stim.freqTrajectory)/1e3+.15], 'LineWidth', lw2);
box off;
set(gca, 'YTick', ytick_vals, 'YTickLabel', ytick_labs);

ax(3)= subplot(223);
spectrogram(y, blackman(anl.tWindow *stim.fs), round(anl.tWindow*anl.tWindow*anl.fracOVlap), anl.nfft, 'centered', 'yaxis', stim.fs);
hold on;
lformant_Han(1)= line(stim.t, 0*stim.t+.15, 'LineWidth', lw);
lformant_Han(2)= line(stim.t, 0*stim.t-.15, 'LineWidth', lw);
set(lformant_Han, 'color', helper.get_color('r'), 'LineStyle', '--');
ylim(yl_val);
set(gca, 'YTick', ytick_vals, 'YTickLabel', ytick_labs);
colorbar off;
ylHan= ylabel('Frequency (kHz)');
ylHan.Position([1 2])= [-.15 3];
xlabel('Time (s)')

ax(4)= subplot(224);
[Py, Fy]= helper.plot_dft(y, stim.fs, 'sided', 2, 'dc', true, 'plot', false);
plot(Py, Fy/1e3, 'LineWidth', lw);
hold on;
plot([0 0], [-.1 .1], 'LineWidth', lw2);
box off;
xlabel('|DFT| in dB');
set(gca, 'YTick', ytick_vals, 'YTickLabel', ytick_labs);


linkaxes(ax, 'y')
ylim(ax(1), yl_val);
linkaxes(ax([2 4]), 'x')
xlim(ax(4), [max([Px(:);Py(:)])-60 10+max([Px(:);Py(:)])])

set(findall(gcf,'-property','FontSize'),'FontSize', 9);
helper.add_subplot_letter(2, 2, 11, .05, 1.05);

%% define new axes for AB
Xshift_X= .08;
Xwidth_X= .425;
Ywidth_X= .37;
Xcorner_X= .06;
Yshift_X= .11;
Ycorner_X= .105;

% C
set(ax(3),'Position',[Xcorner_X Ycorner_X Xwidth_X Ywidth_X])
drawnow
% A
set(ax(1),'Position',[Xcorner_X Ycorner_X+Ywidth_X+Yshift_X Xwidth_X Ywidth_X])
drawnow
% define new axes for AB
Xwidth_CD=Xwidth_X;
Ywidth_CD= Ywidth_X;
Xcorner_CD= Xcorner_X+Xwidth_X+Xshift_X;
Yshift_CD= Yshift_X;
Ycorner_CD= Ycorner_X;
% Yshift=0.06;

% D
set(ax(4),'Position',[Xcorner_CD Ycorner_CD Xwidth_CD Ywidth_CD])
drawnow
% B
set(ax(2),'Position',[Xcorner_CD Ycorner_CD+Ywidth_CD+Yshift_CD Xwidth_CD Ywidth_CD])
drawnow

if saveFig
    saveas(gcf, [LatexDir 'Fig11'], 'epsc');
end