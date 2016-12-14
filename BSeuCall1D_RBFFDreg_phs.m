function [u,err,tim,x,dx,N,W] = BSeuCall1D_RBFFDreg_phs(N,n,ep,M,Kmul)
%% 1D European Call RBF-FD
% 2016-02-06

tic
%% Parameters
K=100;
Kx=1; %strike
T=1; %maturation
r=0.03; %interest
sig=0.15; %volatility

%% Grid
x=transpose(linspace(0,Kmul,N));
dx=x(2)-x(1);

indin=2:N-1;

dt=T/(M-1);
% t=T:-dt:0;

%% Initial condition
u=max(x-Kx,zeros(N,1)); %u0=u;

%% RBF
phi = 'phs';
m = floor(n/4) %number of polynomial terms;

parallel = 1;
W = BSweights1Drbffd_phs(r,sig,x,N,n,indin,phi,ep,m,parallel);

%% Integration
I=speye(N);

%BDF-1
A=I-W*dt;

u1=u;

b=u1;
b(end)=x(end)-Kx*exp(-r*dt);

u=A\b;
u=max(u,0);

%BDF-2
A=I-(2/3)*dt*W;
rcm=symrcm(A);
A=A(rcm,rcm);
[L1, U1]=lu(A);
for ii=3:M
    u2=u1;
    u1=u;

    b=((4/3)*u1-(1/3)*u2);
    b(end)=x(end)-Kx*exp(-r*(ii-1)*dt);

    u(rcm)=L1\b(rcm);
    u(rcm)=U1\u(rcm);
    u=max(u,0);
end
tim=toc;

%% Error
% indreg=[];
% for ii=1:length(x)
%     %         if (xfd(ii)-1)^2/((0.95*K)^2)+(yfd(ii)-1)^2/((0.95*K)^2)<=1
%     if x(ii)>=1/3*Kx && x(ii)<=5/3*Kx
%         indreg=[indreg ii];
%     end
% end
% 
% x=x(indreg);
% u=u(indreg);

% K=100;
x=K*x;
u=K*u; %u0=K*u0;

ua=rsol(sig, r, K, T, x);

% figure()
% plot(x,u0,'k--',x,u,'b-',x,ua,'r-')

err=(u-ua);

% figure()
% plot(x,abs(u-ua));


end
