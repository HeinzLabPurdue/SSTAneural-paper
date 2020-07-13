function x = loadpic(picNum)     % Load picture
picSearchString = sprintf('p%04d*.m', picNum);
picMFile = dir(picSearchString);
if ~isempty(picMFile)
    run(picMFile.name);
    x= ans;
else
    error('Picture file p%04d*.mat not found.', picNum);
    x = [];
    return;
end
