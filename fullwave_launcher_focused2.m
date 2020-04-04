%addpath /celerina/gfp/mfs/2020_ultrasound_seminar/fullwave/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% WRITTEN: JUNE 21, 2018
% LAST MODIFIED: APRIL 4, 2020
% Launch Fullwave 2 code, easy matlab wrapper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load vishuman_abdominal_slice
figure(1), imagesc(cmap'), figure(2), imagesc(rhomap')

%%% Basic variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c0 = 1540;         % speed of sound (m/s)
omega0 = 2*pi*3.7e6; % center radian frequency of transmitted wave
wX = 8e-2;         % width of simulation field (m)
wY = 14e-2;         % depth of simulation field (m)
duration = wY*2.3/c0;  % duration of simulation (s)
p0 = 1e5; % pressure in Pa
%%% Advanced variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ppw = 12;           % number of points per spatial wavelength
cfl = 0.5;         % Courant-Friedrichs-Levi condition
%%% Grid size calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda = c0/omega0*2*pi;
nX = round(wX/lambda*ppw);  % number of lateral elements
nY = round(wY/lambda*ppw);  % number of depth elements
nT = round(duration*c0/lambda*ppw/cfl);
dX = c0/omega0*2*pi/ppw
dY= c0/omega0*2*pi/ppw
dT = dX/c0*cfl;

ncycles = 2; % number of cycles in pulse
dur = 2; % exponential drop-off of envelope
%%% Generate input coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%%
inmap = zeros(nX,nY); 
%inmap(:,1) = ones(nX,1); inmap(:,2) = ones(nX,1); inmap(:,3) = ones(nX,1);
%inmap(round(nX/2),round(nY/2))=1;
rad=4.9e-2;
xducerrad=4.9e-2;
xducercen=[nX/2 -rad/dX+12.4e-3/dX];
idi=circleIdx(size(inmap),xducercen,xducerrad/dY);
inmap(idi)=1;
imagesc(inmap')

for i=1:nX
    j=find(inmap(i,:)==0); j=j(1);
    inmap(i,1:max([j-8 0]))=0;
end
incoords = mapToCoords(inmap);
[tmp idcr]=sort(incoords(:,1));
incoords=incoords(idcr,:);
plot(incoords(:,1),incoords(:,2),'.')

ss=zeros(nX,nY);
cen = [nX/2 -rad/dX+12.4e-3/dX];

theta=atan2(incoords(:,2)-cen(2),incoords(:,1)-cen(1));
Dtheta=max(theta)-min(theta);
for tt=1:128
    idtheta=find(theta<=max(theta)-Dtheta*(tt-1)/128 & theta>max(theta)-Dtheta*(tt)/128);
    plot(incoords(:,1),incoords(:,2),'.'), hold on
    plot(incoords(idtheta,1),incoords(idtheta,2),'r.'); hold off, drawnow
    incoords(idtheta,3)=tt;
end



foc2=10e-2;

focs=0;
Dtheta=30*pi/180/(128-1);
for tt=1:128
    theta=(tt-(128-1)/2)*Dtheta;
    xx=(foc2+xducerrad)/dX*sin(theta)+xducercen(1);
    yy=(foc2+xducerrad)/dX*cos(theta)+xducercen(2);
    focs(tt,1)=xx;
    focs(tt,2)=yy;
end
plot(focs(:,1),focs(:,2),'g.')
axis equal
min(focs(:,1)) % if this number is negative you are focusing outside the sim

%%% Generate output coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%
%outmap=zeros(nX,nY); modX=1; modY=1;
%[modidy modidz]=meshgrid(1:modX:nX,1:modY:nY);
%outmap(modidy,modidz)=1;
%imagesc(outmap');
%outcoords=mapToCoords(outmap);

%outcoords(:,3)=0;

outmap=zeros(nX,nY);
for i=1:nX
    j=find(inmap(i,:)==1); 
    if(~isempty(j))
        j=j(end)+2;
        outmap(i,j)=1;
    end
end
outcoords=mapToCoords(outmap);
outcoords(:,3)=1;
%outcoords=[outcoords; outcoords2;];
[tmp idcr]=sort(outcoords(:,1));
outcoords=outcoords(idcr,:);

plot(incoords(:,1),incoords(:,2),'.'), hold on
plot(outcoords(:,1),outcoords(:,2),'r.')
plot(focs(:,1),focs(:,2),'g.'), hold off
axis equal


%%% Generate field maps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
materials
rho0=1000;
A0 = 3/(20/log(10));
N0 = -tissue.beta/(rho0*c0^4);
fnumber=foc/nY;
beamwidth=round(lambda*fnumber/dY);
nlines=1;
[m.c m.rho m.A m.beta] = img2fieldFlatten2('r102gh.tif',dX,dY);


cmap=interp2easy(cmap,dm/dX,dm/dX,'nearest');
rhomap=interp2easy(rhomap,dm/dX,dm/dX,'nearest');
Amap=interp2easy(Amap,dm/dX,dm/dX,'nearest');
Nmap=interp2easy(Nmap,dm/dX,dm/dX,'nearest');

cmap=cmap(:,round(1.2e-2/dX):round(7e-2/dX));
rhomap=rhomap(:,round(1.2e-2/dX):round(7e-2/dX));
Amap=Amap(:,round(1.2e-2/dX):round(7e-2/dX));
Nmap=Nmap(:,round(1.2e-2/dX):round(7e-2/dX));



nXextend=size(cmap,1);
res_cell = rescell2d(c0,omega0,foc*dY,wY,ncycles,dX,dY);
num_scat = 12;
scat_size = dY;
scat_size_samp = round(scat_size/dY);
cscat_lesion = generate_c_scat(1,1,num_scat/res_cell, scat_size_samp, nXextend, nY)-1;
%idl=circleIdx(size(cscat_lesion),[nYextend/2 foc],5e-3/dY);
%cscat_lesion(idl)=0;
%cscat_lesion(:,1:10)=0;
imagesc(cscat_lesion'), colorbar
csr=0.06; % scatterer impedance contrast 
cscat_lesion;


for i=1:size(cmap,1)
    idcmap=find(cmap(i,:)<=fat.c0);idcmap=idcmap(end);
    cmap(i,idcmap:end)=liver.c0;
    rhomap(i,idcmap:end)=liver.rho0;
    Amap(i,idcmap:end)=liver.alpha;
    Nmap(i,idcmap:end)=0;
    cscat_lesion(i,idcmap:end)=0;
end

%cmap=(cmap-c0)/5+c0;
%Amap=Amap*0+0.5;

gfilts=((1:10)/10).^2*ppw/2;
gfilt=(5/10)^2*ppw/2;
rhofacs=0:0.1:1;
rhofac=0;

orig = [round(2e-2/dX) 1]
c = chopField(cmap,liver.c0,orig,nX,nY);
rho = chopField(rhomap,liver.rho0,orig,nX,nY);
A = chopField(Amap,liver.alpha,orig,nX,nY);
cs=chopField(cscat_lesion,0,orig,nX,nY);
beta=c*0;


c=imgaussfilt(c,gfilt);
rho=imgaussfilt(rho,gfilt);
A=imgaussfilt(A,gfilt);
beta=imgaussfilt(beta,gfilt);


crho=c.*rho; crhohomog=crho*0+mean(mean(crho)); rho2=crhohomog./c;
rho=rho*(1-rhofac)+rho2*rhofac;

c=c-cs*c0*csr;

for i=1:nX
    j=1;
    j=find(inmap(i,:)==1);
    if(~isempty(j))
        j=j(end);
    end
    j=j+round(ppw/2);
    c(i,j:end)=c(i,1:end-j+1); c(i,1:j)=connective.c0;
    rho(i,j:end)=rho(i,1:end-j+1); rho(i,1:j)=connective.rho0;
    A(i,j:end)=A(i,1:end-j+1); 
    beta(i,j:end)=beta(i,1:end-j+1); 
end

foc=round(10e-2/dX)+max(incoords(:,2));
fcen=[round(nX/2) foc]; % center of focus
idl=circleIdx(size(c),fcen,2); c(idl)=c0*0.5; rho(idl)=rho0*0.5;

figure(1), clf, imagesc(c'), figure(2), clf, imagesc(rho')
drawnow
  

for tt=1:128
    tt
    idtt=find(incoords(:,3)==tt);
    
    %%% Generate initial conditions based on input coordinates %%%%%%
    t = (0:nT-1)/nT*duration-ncycles/omega0*2*pi;
    icvec = exp(-(1.05*t*omega0/(ncycles*pi)).^(2*dur)).*sin(t*omega0)*p0;
    plot(icvec)
    %icmat=ones(length(idtt),1)*icvec;
    icmat = focusCoords(focs(tt,1),focs(tt,2),incoords,icvec,cfl);
    imagesc(icmat), drawnow
    
    cwd=pwd; addpath(cwd);
    outdir=['fa/fa_' num2str(tt)]; eval(['!mkdir -p ' outdir]); 
    eval(['!cp fullwave2_try6_nln_relaxing ' outdir]);
    cd(outdir)
    launch_fullwave2_try6_nln_relaxing4(c0,omega0,wX,wY,duration,p0,ppw,cfl,c,rho,A,beta,incoords,outcoords,icmat)
    %eval('!./fullwave2_try6_nln_relaxing &')
    cd(cwd);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LAUNCH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ncoordsout=size(outcoords,1);
 
for tt=1:11
    tt
    cwd=pwd; addpath(cwd);
    outdir=['fa/fa_' num2str(tt)]; eval(['!mkdir -p ' outdir]); 
    cd(outdir)
    eval('!./fullwave2_try6_nln_relaxing &')
    cd(cwd);
end

pause(600);

for tt=12:128    
    tt
    outdir=['fa/fa_' num2str(tt-11)]; 
    nRun=sizeOfFile([outdir '/genout.dat'])/4/ncoordsout
    while (nRun<nT-3)
        nRun=sizeOfFile([outdir '/genout.dat'])/4/ncoordsout
        pause(600);
    end
    cwd=pwd; addpath(cwd);
    outdir=['fa/fa_' num2str(tt)]; eval(['!mkdir -p ' outdir]); 
    cd(outdir)
    eval('!./fullwave2_try6_nln_relaxing &')
    cd(cwd);
    pause(10);
end

tt=22
outdir=['fa/fa_' num2str(tt)]; 
nRun=sizeOfFile([outdir '/genout.dat'])/4/ncoordsout
pxducer=readGenoutSlice([outdir '/genout.dat'],0:nRun-2,size(outcoords,1));
imagesc(powcompress(pxducer,1/4))
 
% save fullwave_launcher_workspace
