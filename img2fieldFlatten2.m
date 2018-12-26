function [c rho A beta] = img2fieldFlatten2(fname,dY,dZ)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% WRITTEN: NOV 13, 2013
% LAST MODIFIED: DECEMBER 18, 2018
% read abdominal image and flatten
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

img = imread(fname,'TIFF');


img1=zeros(size(img))+double(max(max(img)));
for j=1:size(img,2)
  vec=double(img(:,j));
  idv=find(vec<200);
  if(~isempty(idv))
    img1(1:end-idv(1)+1,j)=vec(idv(1):end);
  end
end


yaxisimg = (0:size(img1,1)-1)*1e-3/11.81102;
zaxisimg = (0:size(img1,2)-1)*1e-3/11.81102;
[X Y] = meshgrid(yaxisimg,zaxisimg);
[XI YI] = meshgrid(0:dY:yaxisimg(end),0:dZ:zaxisimg(end));
img2 = interp2(X,Y,double(img1'),XI,YI,'nearest')';

materials

c = zeros(size(img2));
alpha = zeros(size(img2));
rho = zeros(size(img2));
beta = zeros(size(img2));

% water/liver
[row col] = find(img2==250);
for i=1:length(row)
  c(row(i),col(i)) = 1540;
  rho(row(i),col(i)) = 1000;
  alpha(row(i),col(i)) = 3;
  beta(row(i),col(i)) = tissue.beta;
end
% fat
[row col] = find(img2==200);
for i=1:length(row)
  c(row(i),col(i)) = 1478;
  rho(row(i),col(i)) = 950;
  alpha(row(i),col(i)) = 1.8;
  beta(row(i),col(i)) = fat.beta;
end

% muscle
[row col] = find(img2==138);
for i=1:length(row)
  c(row(i),col(i)) = 1547;
  rho(row(i),col(i)) = 1050;
  alpha(row(i),col(i)) = 4.1;
  beta(row(i),col(i)) = muscle.beta;
end

% connective tissue
[row col] = find(img2==0);
for i=1:length(row)
  c(row(i),col(i)) = 1613;
  rho(row(i),col(i)) = 1120;
  alpha(row(i),col(i)) = 5.9;
  beta(row(i),col(i)) = skin.beta;
end


%A = 2*c0^2*(alpha)/omega0^2/c0^4;
A = alpha/(20/log(10));

c=c';
A=A';
rho=rho';
beta=beta';
