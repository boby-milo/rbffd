% function [u,err,tim,x,dx,N,W] = HSTeuCall1D_RBFFDreg(N,n,ep,M)
%K,T,r,sig,M
%% 1D European Call RBF-FD
% Copyright 2015, Slobodan Milovanovic
% 2016-02-06
clc
close all
clear

n=25;

phi='gs';
ep=1;


Nx=300;
Ny=round(Nx/8);

S = [90,100,110];
V = 0.0225;

%U = [2.302535842814927 7.379832496149447 14.974005277144057]; %ref solution
% par = {S,K,T,r,V,kap,th,sig,rho};

tic
%% Parameters
K = 100; %strike
Kx = 1;
T = 1.0;
r = 0.03;

kap = 2;
th = 0.0225;
sig = 0.25;
rho = -0.5;

xmax = 8*Kx;
ymax = 1;

%% Grid
Kx = 1;
x = transpose(linspace(0,xmax,Nx));
y = transpose(linspace(0,ymax,Ny));
dx = x(2)-x(1);
dy = y(2)-y(1);

[X,Y] = meshgrid(x,y);

xvec = transpose(X(:)); yvec = transpose(Y(:));

% ball = xvec+8*yvec;
% ind = find(ball<=8+0.5*dx);

% xvec = xvec(ind);
% yvec = yvec(ind);
% ball = ball(ind);

%indices
ind = 1:numel(xvec);
% indcf = 1;
indcf = find(xvec<=min(xvec)+0.99*dx); indcf = unique(indcf);

%indff1 = find(xvec>=max(xvec)-0.99*dx); indff2 = find(yvec>=max(yvec)-0.99*dy);
%indff = [indff1, indff2]; indff = unique(indff);
indff = find(xvec>=max(xvec)-0.99*dx); indff = unique(indff);

indin = ind; indin([indff,indcf]) = [];

L = 8;

M = 1000;
dt = T/(M-1);
t = T:-dt:0;

N=numel(xvec);

%% Initial condition
u0 = max(xvec-Kx,zeros(1,length(xvec)));
u = u0';

figure(1)
clf
plot(xvec,yvec,'.')
hold on
plot(xvec(indff),yvec(indff),'r^')
plot(xvec(indcf),yvec(indcf),'bs')
plot(xvec(indin),yvec(indin),'k.')
axis equal
axis tight
hold off

figure(2)
tri = delaunay(xvec',yvec');
trisurf(tri, xvec', yvec', u);
shading interp
colorbar
% view(2)
axis vis3d
% axis equal
axis tight


%% RBF
s = [xvec' yvec'];

% Weights
indc = findKNearestNeighbors(s,s,n);

iind = repmat(indin,n,1); iind = iind(:); %n*N
jind = transpose(indc(indin,:)); jind = jind(:);%n*N
Wval = zeros(n,numel(indin));  %n*N

% internal points {
bb = 0;
for ii = indin
    bb = bb+1;
    
    sc = [xvec(ii),yvec(ii)]; xc = sc(:,1); yc = sc(:,2);
    se = s(indc(ii,:),:);
    
    Rc = xcdist(se,se,1);
    A = RBFmat(phi,ep,Rc,'0',1);
    
    Ax = RBFmat(phi,ep,Rc,'1',1);
    Ay = RBFmat(phi,ep,Rc,'1',2);
    
    Axx = RBFmat(phi,ep,Rc,'2',1);
    Ayy = RBFmat(phi,ep,Rc,'2',2);
    Axy = RBFmat(phi,ep,Rc,'m2',1:2);
    
    l = transpose(-r*A(1,:)...
        + r*xc'.*Ax(1,:)...
        + kap*(th-yc').*Ay(1,:)...
        + 0.5*yc'.*xc'.^2.*Axx(1,:)...
        + 0.5*sig^2*yc'.*Ayy(1,:)...
        + rho*sig*xc'.*yc'.*Axy(1,:));
    
    wc = A\l;
    
    Wval(:,bb) = wc;
end
% } internal points
Wval = Wval(:);
W = sparse(iind,jind,Wval,N,N);

%         display('Weights completed');
I=speye(size(W));
%
%% Integration
% BDF-1
u1=u;
A=I-dt*W;

b=u1;
b(indff)=xvec(indff)-Kx*exp(-r*dt);

u=A\b;
u=max(u,zeros(size(u)));

% BDF-2
A=I-(2/3)*dt*W;
rcm=symrcm(A);
A=A(rcm,rcm);
[L1, U1]=lu(A);
for ii=3:M
    u2=u1;
    u1=u;
    b=(4/3)*u1-(1/3)*u2;
    b(indff)=xvec(indff)-Kx*exp(-r*(ii-1)*dt);

    u(rcm)=L1\b(rcm);
    u(rcm)=U1\u(rcm);

    u=max(u,zeros(size(u)));
end
tim=toc;


%% Error
% indreg=[];
% for ii=1:length(xvec)
%     %         if (xfd(ii)-1)^2/((0.95*K)^2)+(yfd(ii)-1)^2/((0.95*K)^2)<=1
%     if xvec(ii)>=1/3*Kx && xvec(ii)<=5/3*Kx && yvec(ii)>=1/3*Kx && yvec(ii)<=5/3*Kx
%         indreg=[indreg ii];
%     end
% end
% 
% xvec=xvec(indreg);
% yvec=yvec(indreg);
% u=u(indreg);
% 
% x=[xvec' yvec'];

xvec=K*xvec;
u=u*K;

figure(3)
tri = delaunay(xvec',yvec');
trisurf(tri, xvec', yvec', u);
shading interp
colorbar
% view(2)
axis vis3d
% axis equal
axis tight



x=[90,100,110];
y=th;
U = [2.302535842814927 7.379832496149447 14.974005277144057];



uinterp=griddata(xvec,yvec,u,x,y,'cubic');

err=max(abs(uinterp-U)./U)


