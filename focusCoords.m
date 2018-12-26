function [icmat] = focusCoords (idy,idz,coords,icvec,cfl)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% WRITTEN: NOV 13, 2013
% LAST MODIFIED: NOV 13, 2013
% Focus coordinates 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dd = sqrt((coords(:,1)-idy).^2+(coords(:,2)-idz).^2);
dd = -(dd/cfl);
dd = dd-min(dd);

icmat = zeros(size(coords,1),length(icvec));
for i=1:size(coords,1)
  icmat(i,dd(i)+1:end)  = icvec(1:end-dd(i));
end
