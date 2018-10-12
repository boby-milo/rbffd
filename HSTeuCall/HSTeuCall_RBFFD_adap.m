function [ures,errres,tim,Acond,N] = HSTeuCall_RBFFD_adap(Nx,pl)
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
yres = ones(size(xres))*0.0225;

dx = 1/Nx;

%% Model
%%% HST
alpha = 0;
beta = 1;
gamma = 0;

fmod = @(x) 0.5*alpha*x.^2 + beta*x + gamma;

%parset
T = 1;
r = 0.03;
kappa = 2;
eta = 0.0225;
sig = 0.25;
rho = -0.5;

K = 100;

%% Grid
Kmul = 4;
Kx = 1/Kmul;

s_max = Kmul*K;
v_max = 0.5;

slim = [0, 1, 0.001, 1];
xscale = s_max/slim(2);
yscale = v_max/slim(4);

sres = [xres'./xscale, yres'./yscale];

[s,N,indres,indin,indcf,indff,inddy0] = HSTeuCall_grid_adap(Nx, Kx, sres, slim, pl);
xvec = s(:,1); yvec = s(:,2);
Nx = sqrt(N)

tri = delaunay(xvec, yvec);
trindin = delaunay(xvec(indin), yvec(indin));

tic
%% Nearest neighbour
indx = findKNearestNeighbors(s,s,n);

par = [-r*ones(size(xvec))...
    r*xvec...
    (1/yscale)*kappa*(eta-yscale*yvec)...
    0.5*beta^2*yscale*xvec.^2.*yvec...
    rho*sig*beta*xvec.*yvec...
    0.5*sig^2*(1/yscale)*yvec];

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
u0 = max(xvec-Kx, 0);
u = u0;

[u,xvec,yvec,Acond] = QLSVeuCallintegrate(u,I,T,M,W,xscale,yscale,xvec,yvec,indin,indff,inddy0,Kx,r,'gmres');
tim = toc;

%% Reference Solution
uref = nan(size(xvec));
for ii = 1:numel(xvec)
    uref(ii) = HSTeuCall_COS(xvec(ii),K,T,r,yvec(ii),kappa,eta,sig,rho);
end

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
% xvec(indres);
% yvec(indres);
% ures = u(indres)';

ures = griddata(xvec,yvec,u,xres,yres,'cubic');

% U(1) = 0.908502728459621;
% U(2) = 9.046650119220969;
% U(3) = 28.514786399298796;
U = [2.302535842814927, 7.379832496149447, 14.974005277144057];

errres = ures-U;
err = u-uref;

end
