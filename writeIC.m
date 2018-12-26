function [] = writeIC (fname,icmat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% WRITTEN: NOV 13, 2013
% LAST MODIFIED: NOV 13, 2013
% write initial condition matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nullvec = zeros(0);
fid = fopen(fname,'wb');
fwrite(fid,nullvec);
fclose(fid);

fid = fopen(fname,'a');
for i=1:size(icmat,2)
  fwrite(fid,icmat(:,i),'float');
end
fclose(fid);
