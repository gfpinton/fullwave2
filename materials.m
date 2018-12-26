%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% WRITTEN: NOV 13, 2013
% LAST MODIFIED: NOV 13, 2013
% common acoustic parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fat = struct('bovera',9.6,'alpha',0.4,'ppower',1.1,'c0',1479,'rho0',937);
fat.beta = 1+fat.bovera/2;
liver = struct('bovera',7.6,'alpha',0.5,'ppower',1.1,'c0',1570,'rho0',1064);
liver.beta = 1+liver.bovera/2;
muscle = struct('bovera',9,'alpha',0.15,'ppower',1.0,'c0',1566,'rho0',1070);
muscle.beta = 1+muscle.bovera/2;
water = struct('bovera',5,'alpha',0.005,'ppower',2.0,'c0',1480,'rho0',1000);
water.beta = 1+water.bovera/2;
skin = struct('bovera',8,'alpha',2.1,'ppower',1,'c0',1498,'rho0',1000);
skin.beta = 1+skin.bovera/2;
tissue = struct('bovera',9,'alpha',0.5,'ppower',1,'c0',1540,'rho0',1000);
tissue.beta = 1+tissue.bovera/2;
connective = struct('bovera',8,'alpha',0.5,'ppower',1,'c0',1613,'rho0',1120);
connective.beta = 1+connective.bovera/2;
blood = struct('bovera',5,'alpha',0.005,'ppower',2.0,'c0',1520,'rho0',1000);
blood.beta = 1+blood.bovera/2;


