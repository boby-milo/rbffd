function [u,err,tim,x,dx,n,N,W] = BSeuCall2Dbasket_RBFFDadap_phs(Nx,p,d,M,Kmul,g)
%% 2D EU Call RBF-FD with BDF2
% 2016-02-04 sparse
load('UrefEU.mat')

tic
%% Model
r=0.03;
sig1=0.15;
sig2=0.15;
rho=0.5; %rho=[1 0.5; 0.5 1];

T=1;
K=1;
%% Grid
Kx=1/Kmul;

% Nx=100;
i=1:Nx;
Ki=2*Kx;
S=1;

% g=5; %tune this! 1,2,3,4,5

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

L=1;
Nlinsq=sqrt(ceil(N*2-sqrt(N*2)));
dx=L/(Nlinsq-1);

% M=100000;
dt=T/(M-1);
t=T:-dt:0;

%% Initial condition
u0=max((1/2)*(xvec+yvec)-Kx,zeros(1,length(xvec)));
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

%% RBF
phi='phs';

dim = 2; %problem dimension

m = nchoosek(p+dim, p); %number of polynomial terms;
n = round(4*m);
s = [xvec' yvec'];

parallel = 0;
W = BSweights2Drbffd_phs(r,sig1,sig2,rho,s,N,n,m,p,indin,phi,d,'reg',parallel);

%% Integration
I = speye(size(W));
% BDF-1
u1 = u;
A = I-dt*W;

b = u1;
b(indff) = 0.5*(xvec(indff)+yvec(indff))-Kx*exp(-r*dt);

u = A\b;
u = max(u,zeros(size(u)));

% BDF-2
A = I-(2/3)*dt*W;
rcm = symrcm(A);
A = A(rcm,rcm);
[L1, U1] = lu(A);
for ii = 3:M
    u2 = u1;
    u1 = u;
    b = (4/3)*u1-(1/3)*u2;
    b(indff) = 0.5*(xvec(indff)+yvec(indff))-Kx*exp(-r*(ii-1)*dt);

    u(rcm) = L1\b(rcm);
    u(rcm) = U1\u(rcm);

    u = max(u,zeros(size(u)));
end
tim = toc;


%% Error
% indreg = [];
% for ii = 1:length(xvec)
%     %         if (xfd(ii)-1)^2/((0.95*K)^2)+(yfd(ii)-1)^2/((0.95*K)^2)<=1
%     if xvec(ii)>=1/3*Kx && xvec(ii)<=5/3*Kx && yvec(ii)>=1/3*Kx && yvec(ii)<=5/3*Kx
%         indreg = [indreg ii];
%     end
% end

xvec = K*Kmul*xvec;
yvec = K*Kmul*yvec;
dx = K*Kmul*dx;
u = K*Kmul*u;


x = [xvec' yvec'];


uinterp = griddata(xulti,yulti,uulti,xvec,yvec,'cubic');

err = uinterp'-u;

end
