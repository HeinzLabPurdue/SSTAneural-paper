% function FigS2_create_kinVowel_spectrum(saveFig, LatexDir)
function FigS2_create_kinVowel_spectrum(saveFig, LatexDir)

if ~exist('saveFig', 'var')
    saveFig= 0;
end
if ~exist('LatexDir', 'var')
    LatexDir= ['figures' filesep];
end

plt.xtick_psd_val= [10 100 1e3 5e3];
plt.xtick_psd_lab= cellfun(@(x) num2str(x), num2cell(plt.xtick_psd_val/1e3), 'uniformoutput', false);
plt.fSize= 9;
plt.arrow.HeadLength= 5;
plt.arrow.HeadWidth= 7;
plt.lw= 1.1;


figSize_cm= [15 5 13.2 7];
figure_prop_name = {'PaperPositionMode','units','Position', 'Renderer'};
figure_prop_val =  { 'auto'            ,'centimeters', figSize_cm, 'painters'};  % [Xcorner Ycorner Xwidth Ywidth]
figure(1);
clf;
set(gcf,figure_prop_name,figure_prop_val);

[sig, fs]= audioread(['stimuli' filesep 's2_vowelA_kin.wav']);
sig= helper.gen_rescale(sig, 70);


[~,~,lHan]= helper.plot_dft(sig, fs, 'yscale', 'dbspl', 'yrange', 70);
xlim([80 5e3]);
box off;
set(lHan, 'color', helper.get_color('b'), 'linew', plt.lw);
xlabel('Frequency (kHz)');
ylabel('DFT-Magnitude (dB)');

set(gca, 'XTick', plt.xtick_psd_val, 'XTickLabel', plt.xtick_psd_lab);
set(findall(gcf,'-property','FontSize'),'FontSize', plt.fSize);

annotation('arrow',[.175, .205], [.85, .85], 'HeadWidth', plt.arrow.HeadWidth, 'HeadLength', plt.arrow.HeadLength);
annotation('arrow',[.525, .495], [.85, .85], 'HeadWidth', plt.arrow.HeadWidth, 'HeadLength', plt.arrow.HeadLength);
annotation('arrow',[.625, .675], [.7, .7], 'HeadWidth', plt.arrow.HeadWidth, 'HeadLength', plt.arrow.HeadLength);
text(.06, .98, 'F_0', 'Units', 'normalized');
text(.48, .98, 'F_1', 'Units', 'normalized');
text(.65, .8, 'F_2', 'Units', 'normalized');
text(.8, .52, 'F_3', 'Units', 'normalized');

plt.tick_len= [.02 .02];
set(findall(gcf,'-property','TickLength'),'TickLength', plt.tick_len, 'units', 'normalized');

set(gca, 'Units', 'normalized', 'Position', [.08 .14 .9 .8]);

if saveFig
   saveas(gcf, [LatexDir 'FigS2'], 'epsc'); 
end