function [cout] = chopField(c,c0,orig,nY,nZ)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% WRITTEN: NOV 13, 2013
% LAST MODIFIED: NOV 13, 2013
% chop an acoustic field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cout = zeros(nY,nZ)+c0;

idy = orig(1):orig(1)+nY-1;
idz = orig(2):orig(2)+nZ-1;

idy2 = 1:size(c,1);
idz2 = 1:size(c,2);

idy = idy(find(idy<=idy2(end) & idy>=idy2(1)));
idz = idz(find(idz<=idz2(end) & idz>=idz2(1)));

for i=1:length(idy)
  for j=1:length(idz)
    cout(i,j) = c(idy(i),idz(j));
  end
end


