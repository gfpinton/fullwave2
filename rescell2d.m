function [rescell] = rescell2d(c0,omega0,dy,ay,ncycles,dY,dZ);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% WRITTEN: NOV 13, 2013
% LAST MODIFIED: NOV 13, 2013
% number of elements in a resolution cell 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

resy = 2*pi*c0/omega0*dy/(ay);
resz = 2*pi*c0/omega0*ncycles/2;
rescell = resy/dY*resz/dZ
