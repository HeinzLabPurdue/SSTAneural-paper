clear;
clc;

all_mFiles= dir('Fig*.m');
allFiles= {all_mFiles.name}';

saveFig =1;
FigDir= ['figures' filesep];
if ~isfolder(FigDir)
    mkdir(FigDir);
end

for fileVar= 1:length(allFiles)
    fprintf('Running %s... \n', allFiles{fileVar})
    eval([allFiles{fileVar}(1:end-2) '(saveFig, FigDir)']);
end

