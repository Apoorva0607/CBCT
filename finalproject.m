
%%
clear;

D=200;

Nu=256;
Nv=256;
Np=100;
uoff=Nu/2;
voff=Nv/2;
R0=100;
du=40/Nu;
dv=du;


Nx=256;
Ny=256;
Nz=128;
xoff=Nx/2;
yoff=Ny/2;
zoff=Nz/2;

dx=20/Nx;
dy=20/Nx;
dz=20/Nz;

al=0.5;

E=sheppdef3D();
phantom  = phantom3D(xoff,yoff,zoff,dx,dy,dz,Nx,Ny,Nz,E);
figure;imagesc(flipud((phantom(:,:,zoff))),[1 1.04]); colormap gray; axis image; axis off
figure;imagesc(flipud(squeeze(phantom(:,yoff,:))),[1 1.04]); colormap gray; axis image; axis off
figure;imagesc(flipud(squeeze(phantom(xoff,:,:))),[1 1.04]); colormap gray; axis image; axis off
tstart=tic;
G  = cbdataApoor(uoff,voff,du,dv,Nu,Nv,Np,R0,D,E);
telapsed1=toc(tstart);

%%
%Backprojection
tstart=tic;
B = cbbpjApoor(G, uoff, voff, du, dv, R0, D, xoff, yoff, zoff, dx, dy, dz, Nx, Ny, Nz);
telapsed2=toc(tstart);
figure;
imagesc(flipud(del2(B(:,:,zoff))),[-0.04 0.001]); colormap gray; axis image; axis off
figure;
imagesc(flipud(del2(squeeze(B(:,yoff,:)))),[-0.04 0.001]); colormap gray; axis image; axis off;
figure;
imagesc(flipud(del2(squeeze(B(xoff,:,:)))),[-0.04 0.001]); colormap gray; axis image; axis off;
%% 
%Filtering
F  = RampApoor(G,uoff,voff,du,dv,al,D);
%%
%Filtered backprojection
FilteredG = cbbpjApoor(F, uoff, voff, du, dv, R0, D, xoff, yoff, zoff, dx, dy, dz, Nx, Ny, Nz);
figure;imagesc(flipud((FilteredG(:,:,zoff))),[1 1.04]); colormap gray; axis image; axis off;
figure;imagesc(flipud(squeeze(FilteredG(:,yoff,:))),[1 1.04]); colormap gray; axis image; axis off;
figure;imagesc(flipud(squeeze(FilteredG(xoff,:,:))),[1 1.04]); colormap gray; axis image; axis off;


