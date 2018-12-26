function [map2] = extendMap(map,nbdy)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% WRITTEN: NOV 13, 2013
% LAST MODIFIED: NOV 13, 2013
% Extend map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

map2 = zeros(size(map,1)+2*nbdy,size(map,2)+2*nbdy);
% center 
map2(nbdy+1:end-nbdy,nbdy+1:end-nbdy) = map;

% edges
for i=1:size(map,1)
  map2(i+nbdy,1:nbdy) = ones(1,nbdy)*map(i,1);
  map2(i+nbdy,end-nbdy+1:end) = ones(1,nbdy)*map(i,end);
end
for j=1:size(map,2)
  map2(1:nbdy,j+nbdy) = ones(nbdy,1)*map(1,j);
  map2(end-nbdy+1:end,j+nbdy) = ones(nbdy,1)*map(end,j);
end

% corners
map2(1:nbdy,1:nbdy) = map(1,1);
map2(end-nbdy+1:end,end-nbdy+1:end) = map(end,end);
map2(end-nbdy+1:end,1:nbdy) = map(end,1);
map2(1:nbdy,end-nbdy+1:end) = map(1,end);
