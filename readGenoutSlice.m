function [genout] = readGenoutslice(fname,nTvec,ncoordsout,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% WRITTEN: NOV 13, 2013
% LAST MODIFIED: NOV 13, 2013
% Read output, in slices of time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% genout = readGenoutSlice([outdir 'genout.dat'],nTvec,ncoordsout);
%idc = 1:round(ncoordsout/100):ncoordsout;
%[genout] = readGenoutSlice(fname,nTvec,ncoordsout,idc);

optargin = size(varargin,2);

if(optargin) % stride
  idc = varargin{1};
  genout = single(zeros(length(nTvec),length(idc)));
  fid = fopen(fname,'rb'); 
  for i=1:length(nTvec)
    t = nTvec(i);
    fseek(fid,(t*ncoordsout)*4,-1);
    tmp = fread(fid,ncoordsout,'float=>float'); 
    genout(i,:) = tmp(idc);
    fprintf(1,'\b\b\b\b\b%0.3f',i/length(nTvec));
    %imagesc(genout),colorbar,drawnow
  end
  fclose(fid);
  fprintf(1,'\n');
else % complete
  genout = single(zeros(length(nTvec),ncoordsout));
  fid = fopen(fname,'rb'); 
  for i=1:length(nTvec)
    t = nTvec(i);
    fseek(fid,(t*ncoordsout)*4,-1);
    genout(i,:) = fread(fid,ncoordsout,'float=>float'); 
    fprintf(1,'\b\b\b\b\b%0.3f',i/length(nTvec));
  end
  fclose(fid);
  fprintf(1,'\n');
end
