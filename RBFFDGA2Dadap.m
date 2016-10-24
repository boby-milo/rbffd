function [u,err,tim,x,dx,N,W] = RBFFDGA2Dadap(Nx,n,ep,M)
%% 2D EU Call RBF-FD with BDF2

load('Uref.mat')

tic
%% Parameters
phi='gs';

%% Model
r=0.03;
sig1=0.15;
sig2=0.15;
rho=0.5;

par={r,sig1,sig2,rho};

T=1;
K=1;

%% Grid
Kx=1;

% Nx=100;
i=1:Nx;
Ki=2*Kx;
S=4*Ki;

g=3; %tune this! 1,2,3,4,5

c=2*Ki/g;

dxi=(1/Nx)*(asinh((S-Ki)/c)-asinh(-Ki/c));
xi=asinh(-Ki/c)+i*dxi;
x=[0, Ki+c*sinh(xi)];
y=zeros(numel(x),1);
Kind=-x(x<=Ki)+Ki;

xvec=[]; yvec=[];
for ii=1:numel(x)
    xl=linspace(0,x(ii),ii);
    yl=linspace(x(ii),0,ii);

    xvec=[xvec,xl];
    yvec=[yvec,yl];
end

N=numel(xvec);
ind=1:N;

indcf=1;
indff=(N-numel(x)+1):N;
indin=ind; indin([indff,indcf])=[];

L=8;
Nlinsq=sqrt(ceil(N*2-sqrt(N*2)));
dx=L/(Nlinsq-1);

% M=100000;
dt=T/(M-1);
t=T:-dt:0;

%% Initial condition
u0=max(0.5*(xvec+yvec)-Kx,zeros(1,length(xvec)));
u=u0';

% figure(1)
% clf
% plot(xvec,yvec,'.')
% hold on
% plot(xvec(indff),yvec(indff),'*')
% plot(xvec(indcf),yvec(indcf),'^')
% plot(xvec(indin),yvec(indin),'o')
% axis equal
% axis tight
% hold off

% figure(2)
% tri = delaunay(xvec',yvec');
% trisurf(tri, xvec', yvec', u);
% shading interp
% colorbar
% view(2)
% axis vis3d
% % axis equal
% axis tight
% pause()

%% RBF

s=[xvec' yvec'];

% Weights
W=sparse(zeros(N));
indc=findKNearestNeighbors(s,s,n);
% internal points {
for ii=indin
    sc=[xvec(indc(ii,:))',yvec(indc(ii,:))'];
    se=[xvec(ii),yvec(ii)];

    wc=rbfga_weights_BS('BS2',ep,sc,se,par);

    W(ii,indc(ii,:))=wc;
end
% } internal points

I=speye(size(W));

%% Integration
% BDF-1
u1=u;
A=I-dt*W;
[L1, U1]=lu(A);

b=u1;
b(indff)=0.5*(xvec(indff)+yvec(indff))-Kx*exp(-r*dt);

u=L1\b;
u=U1\u;
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
    b(indff)=0.5*(xvec(indff)+yvec(indff))-Kx*exp(-r*(ii-1)*dt);

    u=L1\b;
    u=U1\u;

    u(rcm)=L1\b(rcm);
    u(rcm)=U1\u(rcm);

    u=max(u,zeros(size(u)));
end
tim=toc;
% figure(3) %solution
% tri = delaunay(xvec,yvec);
% trisurf(tri, xvec', yvec', u);
% axis vis3d
% axis tight
% xlabel('S_1');
% ylabel('S_2');
% zlabel('V(S_1,S_2)');
% drawnow;

%% Error
indreg=[];
for ii=1:length(xvec)
    %         if (xfd(ii)-1)^2/((0.95*K)^2)+(yfd(ii)-1)^2/((0.95*K)^2)<=1
    if xvec(ii)>=1/3*Kx && xvec(ii)<=5/3*Kx && yvec(ii)>=1/3*Kx && yvec(ii)<=5/3*Kx
        indreg=[indreg ii];
    end
end

xvec=xvec(indreg);
yvec=yvec(indreg);

x=xvec;
u=u(indreg);

uinterp=griddata(xulti,yulti,uulti,xvec,yvec,'cubic');

err=uinterp'-u;

% errorreg=max(abs(uinterp'-u));
% errornormreg=rms(uinterp'-u);

% display([ep,errorreg]);
% figure(4) %error plot
% tri = delaunay(xvec',yvec');
% trisurf(tri, xvec', yvec', abs(uinterp'-u));
% shading interp
% colorbar
% view(2)
% axis vis3d
% % axis equal
% axis tight
% xlabel('S_1');
% ylabel('S_2');
% zlabel('\Delta V(S_1,S_2)');
end
