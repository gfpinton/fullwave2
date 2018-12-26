function [idx]=circleIdx(dims,cen,rad)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% WRITTEN: NOV 13, 2013
% LAST MODIFIED: NOV 13, 2013
% find indices of a 2D circle 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idx=zeros(1,ceil(2*rad+1)^2); cc=1;
if(length(dims)==2)
  for i=-ceil(rad):ceil(rad)
    for j=-ceil(rad):ceil(rad)
      if(sqrt(i^2+j^2)<=rad)
        %idx(cc)=round(i+cen(1))*dims(1)+round(j+cen(2));
        idx(cc)=round(i+cen(1))+round(j+cen(2))*dims(1);
        cc=cc+1;
      end
    end
  end
else
  Disp('ERROR! length(dim)~=2!!')
end
idx=idx(1:cc-1);
