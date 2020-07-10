function picNum=getPicNum(filename)
% function picNum=getPicNum(filename)
% Created: M. Heinz 18Mar2004
% For: CNexps (GE/MH)
%
% Returns Picture Number from NEL filename 

if ~isempty(filename)
   if strcmp(filename(1),'p')
      INDs=2:min(strfind(filename,'_'))-1;
      picNum=str2double(filename(INDs));
      return;
   end
end
picNum=NaN;


return;
