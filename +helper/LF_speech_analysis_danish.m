function PSDout= LF_speech_analysis_danish(curChinID, curTrack, curUnit, dB_SPL, tStart3, tEnd3, fName, saveFigs, dirStruct, doPlot)

loading_Dir= dirStruct.loading_dir;
savingDir= dirStruct.Root_savingDir;
latexDir= dirStruct.latexDir;

curChinDatafName=dir(sprintf('%s*%d*',loading_Dir,curChinID));

fs= 20e3;
[stim, fsOrg]=audioread(['stimuli' filesep 'FLN_Stim_S_P.wav']);
stim= helper.gen_resample(stim, fsOrg, fs);
durStim=length(stim)/fs;

spike_data=load(sprintf('%s%s',loading_Dir,curChinDatafName.name));
spike_data=spike_data.spike_data;


if ~isnan(dB_SPL)
    validINDs=find(([spike_data.track]==curTrack) & ([spike_data.unit]==curUnit) & ([spike_data.SPL]==dB_SPL) ==1);
else
    validINDs=find(([spike_data.track]==curTrack) & ([spike_data.unit]==curUnit) & isnan([spike_data.SPL]) ==1);
end

inds2use= validINDs(strcmp({spike_data(validINDs).noise}, 'SSN'));

if ~isempty(inds2use)
    
    clean_speech_spikes_pos=[];
    clean_speech_spikes_neg=[];
    
    for indVar=1:length(inds2use)
        clean_speech_spikes_pos=[clean_speech_spikes_pos; spike_data(inds2use(indVar)).SpikeTrains{1,1}];
        clean_speech_spikes_neg=[clean_speech_spikes_neg; spike_data(inds2use(indVar)).SpikeTrains{1,2}]; %#ok<*AGROW>
    end
    
    binRes=.1e-3;
    binEdges=0:binRes:durStim;
    uRate_pos=histcounts(cell2mat(clean_speech_spikes_pos), binEdges);
    uRate_neg=histcounts(cell2mat(clean_speech_spikes_neg), binEdges);
    
    tc_data=spike_data(inds2use(indVar)).TC;
    
    PSDout= helper.plot_snap_fft_natural_vowel(uRate_pos, uRate_neg, 1/binRes, stim, fs, tStart3, tEnd3, tc_data, doPlot);
    
    PSDout.CF_Hz= spike_data(inds2use(1)).CF_Hz;
    PSDout.fs= 1/binRes;
    
    
    if ~exist('latex_fName', 'var')
        latex_fName= sprintf('%sLF_clean_speech_Q%d_t%d_u%d_tStart%.0fms_tEnd%.0fms_%.0fdBSPL', latexDir, curChinID, curTrack, curUnit, tStart3*1e3, tEnd3*1e3, dB_SPL);
    end
    if contains(latex_fName, 'LF_clean_speech_Q355_t5_u8_tStart480ms_tEnd580ms_65dBSPL')
        plt.tick_len= [.05 .05];
        set(findall(gcf,'-property','TickLength'),'TickLength', plt.tick_len, 'units', 'normalized');
        if saveFigs
            fprintf('Saved as %s.eps\n', latex_fName);
            saveas(gcf, [latexDir 'Fig9'], 'epsc');
        end
    end
else
    warning('No match for Q%d/t%d/u%d/SPL%d', curChinID, curTrack, curUnit, dB_SPL);
end