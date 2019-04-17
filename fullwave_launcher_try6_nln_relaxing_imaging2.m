%addpath /celerina/gfp/mfs/2018_ultrasound_seminar/fullwave/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% WRITTEN: JUNE 21, 2018
% LAST MODIFIED: DECEMBER 18, 2018
% Launch Fullwave 2 code, easy matlab wrapper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Basic variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c0 = 1540;         % speed of sound (m/s)
omega0 = 2*pi*1e6; % center radian frequency of transmitted wave
wX = 2.4e-2;         % width of simulation field (m)
wY = 6e-2;         % depth of simulation field (m)
duration = wY*2.2/c0;  % duration of simulation (s)
p0 = 1e5; % pressure in Pa
%%% Advanced variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ppw = 8;           % number of points per spatial wavelength
cfl = 0.5;         % Courant-Friedrichs-Levi condition
%%% Grid size calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda = c0/omega0*2*pi;
nX = round(wX/lambda*ppw);  % number of lateral elements
nY = round(wY/lambda*ppw);  % number of depth elements
nT = round(duration*c0/lambda*ppw/cfl);
dX = c0/omega0*2*pi/ppw
dY = c0/omega0*2*pi/ppw
dT = dX/c0*cfl;
%%% Generate input coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%%
inmap = zeros(nX,nY); 
inmap(:,1) = ones(nX,1); inmap(:,2) = ones(nX,1); inmap(:,3) = ones(nX,1);
%inmap(round(nX/2),round(nY/2))=1;
incoords = mapToCoords(inmap);
%%% Generate initial conditions based on input coordinates %%%%%%
foc=round(4.5e-2/dX);
fcen=[round(nX/2) foc]; % center of focus

ncycles = 2; % number of cycles in pulse
dur = 2; % exponential drop-off of envelope
t = (0:nT-1)/nT*duration-ncycles/omega0*2*pi;
icmat=zeros(size(incoords,1),nT);
icvec = exp(-(1.05*t*omega0/(ncycles*pi)).^(2*dur)).*sin(t*omega0)*p0;
plot(icvec), hold all
icmat(1:size(incoords,1)/3,:) = focusCoords(fcen(1),fcen(2),incoords(1:size(incoords,1)/3,:),icvec,cfl);
t=t-dT/cfl;
icvec = exp(-(1.05*t*omega0/(ncycles*pi)).^(2*dur)).*sin(t*omega0)*p0;
plot(icvec)
icmat(size(incoords,1)/3+1:size(incoords,1)/3*2,:)=focusCoords(fcen(1),fcen(2),incoords(size(incoords,1)/3+1:size(incoords,1)/3*2,:),icvec,cfl);
t=t-dT/cfl;
icvec = exp(-(1.05*t*omega0/(ncycles*pi)).^(2*dur)).*sin(t*omega0)*p0;
plot(icvec), hold off
icmat(size(incoords,1)/3*2+1:size(incoords,1),:)=focusCoords(fcen(1),fcen(2),incoords(size(incoords,1)/3*2+1:size(incoords,1),:),icvec,cfl);
imagesc(icmat)
%%% Generate output coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%
outmap=zeros(nX,nY); modX=1; modY=1;
[modidy modidz]=meshgrid(1:modX:nX,1:modY:nY);
outmap(modidy,modidz)=1;
imagesc(outmap');
outcoords=mapToCoords(outmap);
%%% Generate field maps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
materials
rho0=1000;
A0 = 3/(20/log(10));
N0 = -tissue.beta/(rho0*c0^4);
fnumber=foc/nY;
beamwidth=round(lambda*fnumber/dY);
nlines=21;
[m.c m.rho m.A m.beta] = img2fieldFlatten2('r102gh.tif',dX,dY);
nYextend=size(m.c,1);
res_cell = rescell2d(c0,omega0,foc*dY,wY,ncycles,dX,dY);
num_scat = 12;
scat_size = dY;
scat_size_samp = round(scat_size/dY);
cscat_lesion = generate_c_scat(1,1,num_scat/res_cell, scat_size_samp, nYextend, nY)-1;
idl=circleIdx(size(cscat_lesion),[nYextend/2 foc],5e-3/dY);
cscat_lesion(idl)=0;
cscat_lesion(:,1:10)=0;
imagesc(cscat_lesion'), colorbar
csr=0.02; % scatterer impedance contrast 


n=1
for n=1:nlines
  orig = [round(1.5e-2/dY-(n-(nlines+1)/2)*beamwidth/4) 1]
  c = chopField(m.c,c0,orig,nX,nY);
  rho = chopField(m.rho,rho0,orig,nX,nY);
  A = chopField(m.A,A0,orig,nX,nY);
  beta = chopField(m.beta,5.5,orig,nX,nY); 
  cs = cscat_lesion(round(nYextend/2-(n-(nlines+1)/2)*beamwidth/4-nX/2):round(nYextend/2-(n-(nlines+1)/2)*beamwidth/4-nX/2)+nX-1,:);

  
  c=imgaussfilt(c,1*ppw/12);
  rho=imgaussfilt(rho,1*ppw/12);
  A=imgaussfilt(A,1*ppw/12);
  beta=imgaussfilt(beta,1*ppw/12);
  
  c=c-cs*c0*csr;
  imagesc(c'), drawnow

  cwd=pwd; addpath(cwd);
  outdir=['lesion' num2str(n)]; eval(['!mkdir -p ' outdir]); 
  eval(['!cp fullwave2_try6_nln_relaxing ' outdir]);
  cd(outdir)
  launch_fullwave2_try6_nln_relaxing4(c0,omega0,wX,wY,duration,p0,ppw,cfl,c,rho,A,beta,incoords,outcoords,icmat)
  if(n<nlines)
   eval('!./fullwave2_try6_nln_relaxing &')
  elseif(n==nlines)
  eval('!./fullwave2_try6_nln_relaxing')
  end
 
 cd(cwd);
end

%%% GENERATE IMAGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
deps = 10e-3:0.125e-3/4:foc*dX*1.15;
lats = 0;
fnumber=1;
idc=find(outcoords(:,2)==4);
xducercoords = outcoords(idc,:);

bm=zeros(length(lats),length(deps),nlines);
idps=cell(length(lats),length(deps));

n=round(nlines/2);
outdir=['lesion' num2str(n) '/']
ncoordsout=size(outcoords,1);
nRun=sizeOfFile([outdir 'genout.dat'])/4/ncoordsout;
pxducer=readGenoutSlice([outdir 'genout.dat'],0:nRun-1,size(outcoords,1),idc);
imagesc(powcompress(pxducer,1/4))
px=pxducer(:,round(size(pxducer,2)/2));
[val idt0]=max(abs(hilbert(px)))


for n=1:nlines
  outdir=['lesion' num2str(n) '/']
  ncoordsout=size(outcoords,1);
  nRun=sizeOfFile([outdir 'genout.dat'])/4/ncoordsout;
  while(nRun<nT-1)
    pause(0.1)
  nRun=sizeOfFile([outdir 'genout.dat'])/4/ncoordsout;
  end
  pxducer=readGenoutSlice([outdir 'genout.dat'],0:nRun-1,size(outcoords,1),idc);
  imagesc(powcompress(pxducer,1/4))
  px=pxducer(:,round(size(pxducer,2)/2));
  %[val idt0]=max(abs(hilbert(px)))

  if(n==1)
    idps=cell(length(lats),length(deps));
    for ii=1:length(lats)
      lat=lats(ii);
      for jj=1:length(deps)
        dep=deps(jj);
        fcen=round([lat/dX+mean(xducercoords(:,1)) dep/dX ]);
        idx=find(abs(xducercoords(:,1)-fcen(1))<=fcen(2)/fnumber);
        dd=focusProfile(fcen,xducercoords(idx,:),dT/dY*c0);
        idt=idt0+round(2*dep/double(mean(mean(c)))/(dT));
        idp=double((size(pxducer,1)*(idx-1))+double(idt)+dd);
        idps{ii,jj}=idp;
      end
    end
  end

  for ii=1:length(lats)
    for jj=1:length(deps)
      bm(ii,jj,n)=sum(pxducer(idps{ii,jj}));
    end
  end
end



%% PLOT THE BMODE IMAGE %%
figure(1)
n=1:nlines; bws=((n-(nlines+1)/2)*beamwidth/4)*dY;
imagesc(bws*1e3,deps*1e3,dbzero(abs(hilbert(squeeze(bm)))),[-40 0])
colormap gray
xlabel('mm'), ylabel('mm')
axis equal, axis tight

figure(2)
n=round(nlines/2);
orig = [round(1.5e-2/dY-(n-(nlines+1)/2)*beamwidth/4) 1]
c = chopField(m.c,c0,orig,nX,nY);
cs = cscat_lesion(round(nYextend/2-(n-(nlines+1)/2)*beamwidth/4-nX/2):round(nYextend/2-(n-(nlines+1)/2)*beamwidth/4-nX/2)+nX-1,:);
c=c-cs*c0*csr;

imagesc(((1:nX)-nX/2)*dX*1e3,(1:nY)*dY*1e3,c')
axis equal, axis tight
axis([bws(1) bws(end) deps(1) deps(end)]*1e3)
xlabel('mm'), ylabel('mm')

