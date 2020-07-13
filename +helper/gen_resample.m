% function vecOut=gen_resample(vecIn, fsOld, fsNew)
function vecOut=gen_resample(vecIn, fsOld, fsNew)

%% find p and q for resampling
temp= rats(round(fsNew)/round(fsOld));
temp(temp==' ')= [];
if ~isempty(find(temp=='/', 1))
    pqEndPartition= find(temp=='/');
    p= str2double(temp(1:pqEndPartition-1));
    q= str2double(temp(pqEndPartition+1:end));
else
    p= str2double(temp);
    q= 1;
end

vecOut=resample(vecIn, p, q);