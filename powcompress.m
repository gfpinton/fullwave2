function [img] = powcompress(img,pow)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% WRITTEN: NOV 13, 2013
% LAST MODIFIED: NOV 13, 2013
% Power compression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:size(img,1)
  for j=1:size(img,2)
    if(img(i,j)>0)
      img(i,j) = img(i,j)^pow;
    else
      img(i,j) = -((-img(i,j))^pow);
    end
  end
end

