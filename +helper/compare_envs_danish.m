function compare_envs_danish(curChinID, curTrack, curUnit, dB_SPL, tStart, tEnd, saveFigs, dirStruct)

loading_Dir=dirStruct.loading_Dir;
latexDir= dirStruct.latexDir;

binRes=.5e-3;

curChinDatafName=dir(sprintf('%s*%d*',loading_Dir, curChinID));

[stim, fs]=audioread(['stimuli' filesep 'FLN_Stim_S_P.wav']);
tStim=(1:length(stim))/fs;
inds2use_stim= tStim>tStart & tStim<tEnd;
stimSnap= stim(inds2use_stim);
tStimSnap= tStim(inds2use_stim);
durStim=length(stim)/fs;

spike_data=load(sprintf('%s%s',loading_Dir,curChinDatafName.name));
spike_data=spike_data.spike_data;

if ~isnan(dB_SPL)
    validINDs=find(([spike_data.track]==curTrack) & ([spike_data.unit]==curUnit) & ([spike_data.SPL]==dB_SPL) ==1);
else
    validINDs=find(([spike_data.track]==curTrack) & ([spike_data.unit]==curUnit) & isnan([spike_data.SPL]) ==1);
end

inds2use_stim= validINDs(strcmp({spike_data(validINDs).noise}, 'SSN'));

if ~isempty(inds2use_stim)
    
    clean_speech_spikes_pos=[];
    clean_speech_spikes_neg=[];
    
    for indVar=1:length(inds2use_stim)
        clean_speech_spikes_pos=[clean_speech_spikes_pos; spike_data(inds2use_stim(indVar)).SpikeTrains{1,1}];
        clean_speech_spikes_neg=[clean_speech_spikes_neg; spike_data(inds2use_stim(indVar)).SpikeTrains{1,2}]; %#ok<*AGROW>
    end
    
    binEdges=0:binRes:durStim;
    tCenters=mean([binEdges(1:end-1); binEdges(2:end)], 1);
    
    inds2use_uRate= tCenters>tStart & tCenters<tEnd;
    
    uRate_pos= histcounts(cell2mat(clean_speech_spikes_pos), binEdges);
    uRate_neg= histcounts(cell2mat(clean_speech_spikes_neg), binEdges);
    
    
    uRate_sum= .5*(uRate_pos+uRate_neg);
    uRate_tfs= .5*(uRate_pos-uRate_neg);
    uRate_tfs(isnan(uRate_tfs))= 0;
    uRate_hil= abs(hilbert(uRate_tfs));
    
    % 
    tMax= 860e-3;
    uRate_pos(tCenters>tMax)= nan;
    uRate_neg(tCenters>tMax)= nan;
    uRate_sum(tCenters>tMax)= nan;
    uRate_tfs(tCenters>tMax)= nan;
    uRate_hil(tCenters>tMax)= nan;
    stim(tStim>tMax)= nan;


    Amax= 1.05*max(abs([uRate_pos(inds2use_uRate), uRate_neg(inds2use_uRate)]));
    yMinMax=[-1.6*Amax 1.6*Amax];
    
    lw=1;
    nSProws= 5;
    nSPcols= 1;
    
    figSize_cm= [15 5 13.2 8];
    figure_prop_name = {'PaperPositionMode', 'units', 'Position', 'Renderer'};
    figure_prop_val =  { 'auto', 'centimeters', figSize_cm, 'painters'};  % [Xcorner Ycorner Xwidth Ywidth]
    figure(1);
    clf;
    set(gcf,figure_prop_name, figure_prop_val);
    
    colGray= 2*helper.get_color('gray');
    
    ax(1)=subplot(nSProws, nSPcols, 1:2);
    hold on;
    
    tCenters_ms= tCenters*1e3;
    plot(tCenters_ms, uRate_pos, 'color', colGray);
    plot(tCenters_ms, -uRate_neg, 'color', colGray);
    nl(1)=plot(tCenters_ms(inds2use_uRate), uRate_pos(inds2use_uRate), 'color', helper.get_color('b'), 'linew', lw);
    nl(2)=plot(tCenters_ms(inds2use_uRate), -uRate_neg(inds2use_uRate), 'color', helper.get_color('r'), 'linew', lw);
    
    stimRatio= .3;
    plot(tStim*1e3, -(1+stimRatio)*Amax+stimRatio*Amax*stim/max(stimSnap), 'color', colGray);
    plot(tStimSnap*1e3, -(1+stimRatio)*Amax+stimRatio*Amax*stimSnap/max(stimSnap), 'color', helper.get_color('k'));
    nl(3)= plot(nan, nan, 'color', helper.get_color('k'), 'linew', lw);
    
    txtHan(1)= text(.02, 1.01, 'A. p(t) & n(t)', 'units', 'normalized');
    [lgHan, icons]= legend(nl, {'p(t)', '-n(t)', 'Stim'}, 'Box', 'off');
    
    lgHan.Box= 'off';
    
    icons(4).XData= mean(icons(4).XData) + [0.1 .3];
    icons(4).LineWidth= 1;
    icons(6).XData= mean(icons(6).XData) + [0.1 .3];
    icons(6).LineWidth= 1;
    icons(8).XData= mean(icons(8).XData) + [0.1 .3];
    icons(8).LineWidth= 1;
    
    ylim(yMinMax);
    
    
    set(gca, 'xticklabel', '');
    axis tight;
    
    
    ax(2)=subplot(nSProws, nSPcols, 3);
    hold on;
    plot(tCenters_ms, uRate_sum, 'color', colGray);
    plot(tCenters_ms(inds2use_uRate), uRate_sum(inds2use_uRate), 'color', helper.get_color('b'), 'linew', lw);
    txtHan(2)= text(.02, .95, 'B. s(t)', 'units', 'normalized');
    set(gca, 'xticklabel', '');
    axis tight;
    
    
    ax(3)=subplot(nSProws, nSPcols, 4);
    hold on;
    plot(tCenters_ms, uRate_tfs, 'color', colGray);
    plot(tCenters_ms(inds2use_uRate), uRate_tfs(inds2use_uRate), 'color', helper.get_color('b'), 'linew', lw);
    txtHan(3)= text(.02, .95, 'C. d(t)', 'units', 'normalized');
    set(gca, 'xticklabel', '');
    axis tight;
    
    
    ax(4)=subplot(nSProws, nSPcols, 5);
    hold on;
    plot(tCenters_ms, uRate_hil, 'color', colGray);
    plot(tCenters_ms(inds2use_uRate), uRate_hil(inds2use_uRate), 'color', helper.get_color('b'), 'linew', lw);
    
    txtHan(4)= text(.02, .95, 'D. e(t)', 'units', 'normalized');
    axis tight;
    ylHan= ylabel('Discharge rate (spikes/bin)');
    
    linkaxes(ax, 'x');
    linkaxes(ax(2:4), 'y');
    xlabel('Time (ms)');
    ylim(ax(1), [-1+min(ylim(ax(1))) 2+max(ylim(ax(1)))]);
    ylim(ax(end), [-.5+min(uRate_tfs(inds2use_uRate)) .5+max(uRate_sum(inds2use_uRate))]);
    xlim(ax, [tStart-.015 tEnd+.03]*1e3);
    
    ylHan.Position(1:2)= [tStart*1e3-23 2.7*Amax];
    
    
    set(findall(gcf,'-property','FontSize'),'FontSize', 9);
    set(txtHan, 'FontSize', 11);
    
    %%
    Xcorner_X= .075;
    Xwidth_X= .89;
    Ywidth_X= .15;
    Yshift_X= .03;
    Ycorner_X= .12;
    
    % D
    set(ax(4),'Position',[Xcorner_X, Ycorner_X, Xwidth_X, Ywidth_X])
    drawnow
    % C
    set(ax(3),'Position',[Xcorner_X, Ycorner_X+Ywidth_X+Yshift_X, Xwidth_X, Ywidth_X])
    drawnow
    % B
    set(ax(2),'Position',[Xcorner_X, Ycorner_X+2*Ywidth_X+2*Yshift_X, Xwidth_X, Ywidth_X])
    drawnow
    % A
    set(ax(1),'Position',[Xcorner_X, Ycorner_X+3*Ywidth_X+3*Yshift_X, Xwidth_X, 2*Ywidth_X])
    drawnow
    
    if saveFigs
        if ~exist('latex_fName', 'var')
            latex_fName= sprintf('%senv_comp_Q%d_t%d_u%d_tStart%.0fms_tEnd%.0fms_%.0fdBSPL', latexDir,curChinID, curTrack, curUnit, tStart*1e3, tEnd*1e3, dB_SPL);
        end
        if contains(latex_fName, 'Q355_t2_u4_tStart740ms_tEnd850ms_65dBSPL')
%             fprintf('Saved as %s.eps\n', latex_fName);
            %             saveas(gcf, latex_fName, 'epsc');
            saveas(gcf, [latexDir 'Fig10'], 'epsc');
        end
    end
else
    warning('No match for Q%d/t%d/u%d/SPL%d', curChinID, curTrack, curUnit, dB_SPL);
end