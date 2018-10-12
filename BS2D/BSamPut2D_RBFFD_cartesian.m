function [ures,errres,tim,Acond,N] = BSamPut2D_RBFFD_cartesian(Nx,pl)
if nargin == 1
    pl = 0;
end

%% Parameters
M = 100;

phi = 'phs';
ep = 5;
p = 4;
nm = 5;
dim = 2; %problem dimension
m = nchoosek(p+dim, p) %number of polynomial terms;
n = round(nm*m); n = max(9,n)

xres = [90,100,110];
yres = xres;
dx = 1/Nx;

%% Model
r = 0.03;
sig1 = 0.15;
sig2 = 0.15;
rho = 0.5; %rho=[1 0.5; 0.5 1];

T = 1;
K = 100;

%% Grid
Kmul = 8;
Kx = 1/Kmul;

x_max = Kmul*K;
y_max = x_max;

slim = [0, 1, 0, 1];
xscale = x_max/slim(2);
yscale = y_max/slim(4);

sres = [xres'./xscale, yres'./yscale];

[s,N,indres,indin,indcf,indff] = BSeuCall2D_grid_cartesian(Nx, Kx, sres, slim, pl);
xvec = s(:,1); yvec = s(:,2);
Nx = sqrt(N)

tri = delaunay(xvec, yvec);
trindin = delaunay(xvec(indin), yvec(indin));

tic
%% Nearest neighbour
indx = findKNearestNeighbors(s,s,n);

par = [-r*ones(size(xvec)) ...
    r*xvec ...
    r*yvec ...
    0.5*sig1^2*xvec.^2 ...
    rho*sig1*sig2*xvec.*yvec ...
    0.5*sig2^2*yvec.^2];

[W,wcond] = RBFFDweights(par,s,N,phi,ep,p,n,indin,indx,dim);
I = speye(size(W));
[wcondmax, icondmax] = max(wcond);
indcondmax = find(wcond>=1e16);

if pl
    figure(1)
    plot(xvec(indin(icondmax)),yvec(indin(icondmax)),'*')
    hold on
    plot(xvec(indin(indcondmax)),yvec(indin(indcondmax)),'o')
    drawnow
    
    figure(2)
    trisurf(trindin, xvec(indin), yvec(indin), log10(wcond));
    shading interp
    colorbar
%     view(3)
    axis vis3d
    % axis equal
    axis tight
    view(2)
    drawnow
    
end

%% Integration
u0 = max(Kx-0.5*(xvec+yvec), 0);
u = u0;
lambda=zeros(N,1);
[u,xvec,yvec,Acond] = BSamPut2D_integrate(u,lambda,I,T,M,N,W,xscale,yscale,xvec,yvec,indin,indff,Kx,r,'gmres');

tim = toc;

%% Reference Solution
load('UrefAM.mat')
xulti = 100*xulti;
yulti = 100*yulti;
uulti = 100*uulti;
uref = griddata(xulti,yulti,uulti,xvec,yvec,'cubic');

if pl
    figure(3)
    trisurf(tri, xvec, yvec, u);
    shading interp
    colorbar
    view(3)
    axis vis3d
    % axis equal
    axis tight
    drawnow
    
    figure(4)
    trisurf(tri, xvec, yvec, W*u);
    shading interp
    colorbar
    view(3)
    axis vis3d
    % axis equal
    axis tight
    drawnow
    
    figure(5)
    trisurf(tri, xvec, yvec, log10(abs(u-uref)));
    shading interp
    colorbar
    view(2)
    axis vis3d
    caxis([-6 1])
    axis tight
    drawnow
end

%% Result
% xres = xvec(indres);
% yres = yvec(indres);
ures = griddata(xvec,yvec,u,xres,yres,'cubic');

U = griddata(xulti,yulti,uulti,xres,yres,'cubic');

errres = ures-U;
err = u-uref;

end
