function [] = writeVabs (typ,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% WRITTEN: NOV 13, 2013
% LAST MODIFIED: NOV 13, 2013
% write variables of different types
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% writeVabs('float',dX,'dX',dT,'dT',rho0,'rho0',mu0,'mu0',lambda0,'lambda0',wght,'wght');
% writeVabs('int',nX,'nX',nY,'nY',nT,'nT',nTic,'nTic',modX,'modX',modY,'modY',modT,'modT',npml,'npml',nbdw,'nbdw',ncoordsin,'ncoordsin',ncoordsout,'ncoordsout');
%writeVabs('char',outdir,'outdir');

optargin = size(varargin,2);

% check even number of arguments
if(mod(optargin,2))
  disp('ERROR, optargin not even');
  return;
else
  for k=1:2:optargin
    fname = [varargin{k+1} '.dat'];
    disp(['Writing ' fname])
    fid = fopen(fname,'wb');
    fwrite(fid,varargin{k},typ);
    fclose(fid);
  end
end

