function [imgout] = interp2easy(img,interpx,interpy,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% WRITTEN: NOV 13, 2013
% LAST MODIFIED: NOV 13, 2013
% [imgout] = interp2easy(img,interpx,interpy)
% interpx = 3;
% interpy = 5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
interpstr='*spline';
optargin = size(varargin,2);
if(optargin==1)
  interpstr=varargin{1};
end

[x y] = meshgrid(0:size(img,2)-1,0:size(img,1)-1);
%[xi yi] = meshgrid((0:ceil((size(img,2)*interpy-1)))/interpy,(0:ceil((size(img,1)*interpx-1)))/interpx);
dxi = 1/round((size(img,2))*interpx-1);
xvec = 0:dxi:1;
xvec = xvec*(size(img,2)-1);
dyi = 1/round((size(img,1))*interpy-1);
yvec = 0:dyi:1; 
yvec = yvec*(size(img,1)-1);
[xi yi] = meshgrid(xvec,yvec);

imgout = interp2(x,y,img,xi,yi,interpstr);
%imgout = interp2(x,y,img,xi,yi,'*cubic');
%imgout = interp2(x,y,img,xi,yi,'linear*');

[idi idj] = find(imgout==NaN);
if(isempty(idi))
  %disp('No NaNs found');
else
  disp('Replacing NaNs with zeros');
  imgout(idi,idj) = 0;
end


idi=find(imgout==Inf); imgout(idi)=0;
idi=find(imgout==-Inf); imgout(idi)=0;

