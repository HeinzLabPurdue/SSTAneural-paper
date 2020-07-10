function txtHan= add_subplot_letter(nSProws, nSPcols, fSize, xShift, yShift)

if ~exist('xShift', 'var')
    xShift= 0.01;
end
if ~exist('yShift', 'var')
    yShift= 1.05;
end

SP_letters= 'ABCDEFGHI';
count= 0;
txtHan= nan(nSProws*nSPcols,1);
for rowVar=1:nSProws
    for colVar=1:nSPcols
        count= count+1;
        subplot(nSProws, nSPcols, count);
        txtHan(count)= text(xShift, yShift, ['\bf' SP_letters(count) '\rm'], 'FontSize', fSize, 'Units', 'normalized');
    end
end