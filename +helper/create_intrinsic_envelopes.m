function create_intrinsic_envelopes(curSpikes_pos, curSpikes_neg, anl, filt_obj_LP, filt_abj_BP, fName, dirStruct, saveFig)

[sigEng, fsEng_Org]= audioread(['stimuli' filesep 'Stim_S_P.wav']);

sigEng= helper.gen_resample(sigEng, fsEng_Org, anl.stimFs);
tEng= (1:length(sigEng))/anl.stimFs;

stimDur= length(sigEng)/anl.stimFs;

tEdges= 0:1/anl.fs:stimDur;
tBinCenters= tEdges(1:end-1)+1/2/anl.fs;

uRatePos= histcounts(curSpikes_pos, tEdges);
uRateNeg= histcounts(curSpikes_neg, tEdges);

uRateComp= (uRatePos-uRateNeg)/2;

env_bp= cell(length(anl.BP.freqs), 1);
for fmVar=1:length(anl.BP.freqs)
    cur_filt= filt_abj_BP{fmVar};
    env_bp{fmVar}= filtfilt(cur_filt, uRatePos)';
end

%% plot
Astim= 0.5*max(abs(uRateComp));
Aenv= 1.25*max(abs(cell2mat(env_bp(:))));

figSize_cm= [15 5 13.2 6.5]; % [Xcorner Ycorner Xwidth Ywidth]
fig_num= 1;
figure_prop_name = {'PaperPositionMode','units','Position'};
figure_prop_val =  { 'auto'            ,'centimeters', figSize_cm};  % [Xcorner Ycorner Xwidth Ywidth]
figure(fig_num);
clf;
set(fig_num, figure_prop_name, figure_prop_val);
plt.tick_len= [.025 .025];


sp_ax(1)= subplot(121);
hold on;
lHan(1)= plot(tBinCenters, uRatePos, 'Color', helper.get_color('b'));
lHan(2)= plot(tEng, -1*Astim + Astim*sigEng/max(abs(sigEng)), 'Color', helper.get_color('gray'));


box off;
axis tight;
ylabel('Discharge rate (spikes/bin)');


xlabHan= xlabel('Time (s)');
ylim([-2.5*Astim 1.5*max(uRatePos)]);

sp_ax(2)= subplot(122);
hold on;
txHan= nan(length(anl.BP.freqs), 1);
text(-.37, Aenv*(length(anl.BP.freqs)+1)*1.03, 'F_m (Hz)');
for fmVar=1:length(anl.BP.freqs)
    plot(tBinCenters, env_bp{fmVar} + Aenv*fmVar, 'Color', helper.get_color('b'), 'linew', 1);
    txHan(fmVar)= text(-.27, Aenv*fmVar*1.05, sprintf('%.0f', anl.BP.freqs(fmVar)), 'HorizontalAlignment', 'center');
end
set(gca, 'YTickLabel', '');
xlim([-.55 stimDur]);
ylim([0.5 length(anl.BP.freqs)+1.75]*Aenv);
box off;

[lg, icons]= legend(lHan, 'p(t)', 'Stim');
lg.Box= 'off';lg.FontSize= 18;
hl = findobj(icons,'type','line');
set(hl([1 3]),'LineWidth',1.5);
icons(3).XData= mean(icons(3).XData) + [0 +.25];
icons(5).XData= mean(icons(5).XData) + [0 +.25];


arh=annotation('arrow', [0.44, 0.6],[0.61,0.61]);
arh.LineWidth= 5;
arh.Color= helper.get_color('b');
arh.HeadWidth= 25;

txt= text(-0.49, 0.47, {'Modulation';'Filterbank'}, 'units', 'normalized');
txt.FontWeight= 'bold';
txt.FontSize= 9;
set(findall(gcf,'-property','FontSize'),'FontSize', 9);
helper.add_subplot_letter(1, 2, 11, .05, 1.05);
set(findall(gcf,'-property','TickLength'),'TickLength', plt.tick_len, 'units', 'normalized');
xlabHan.Position(1:2)= [4 -6.7];

% Set axes placement/size
Xwidth=.36;
Xcorner=.072;
Xshift=.18;
Ywidth=.82;
Ycorner=.11;

% A
set(sp_ax(1),'Position',[Xcorner Ycorner Xwidth Ywidth])
drawnow
% B
set(sp_ax(2),'Position',[Xcorner+Xwidth+Xshift Ycorner Xwidth Ywidth])
drawnow

if saveFig && contains(fName, 'HilbEnv_clean_speech_Q277_t2_u1_61dBSPL')
    saveas(gcf, [dirStruct.latexDir 'Fig5'], 'epsc');
end