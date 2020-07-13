function txtHan= add_subplot_letter_txp(nSProws, nSPcols, fSize, xShift, yShift)

if ~exist('xShift', 'var')
    xShift= 0.01;
end
if ~exist('yShift', 'var')
    yShift= 1.05;
end

SP_letters= 'ABCDEFGHI';
txtHan= nan(nSProws*nSPcols,1);
for colVar=1:nSPcols
    for rowVar=1:nSProws
        count_SP_index= colVar + nSPcols*(rowVar-1);
        count_SP_letter= rowVar + nSProws*(colVar-1);
        subplot(nSProws, nSPcols, count_SP_index);
        txtHan(count_SP_index)= text(xShift, yShift, ['\bf' SP_letters(count_SP_letter)], 'FontSize', fSize, 'Units', 'normalized');
    end
end