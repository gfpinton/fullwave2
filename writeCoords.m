function [] = writeCoords(fname,coords)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% WRITTEN: NOV 13, 2013
% LAST MODIFIED: NOV 13, 2013
% write coordinate vectors to file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nullvec = zeros(0);
fid = fopen(fname,'wb');
fwrite(fid,nullvec);
fclose(fid);

fid = fopen(fname,'a');
for i=1:size(coords,2)
  fwrite(fid,coords(:,i),'int');
end
fclose(fid);
