function [u,err,tim,x,dx,N,W] = RBFFDGA1DregS0(N,n,ep,M)
%% 1D European Call RBF-FD
% Copyright 2015, Slobodan Milovanovic

tic
%% Parameters
K=100;
Kx=1; %strike
T=1; %maturation
r=0.03; %interest
sig=0.15; %volatility

%% Grid
% N=8001;
x=transpose(linspace(0,4,N));
dx=x(2)-x(1);

% n=3; %stencil size
m=round((n-1)/2);

% M=100000;
dt=T/(M-1);
t=T:-dt:0;

%% Initial condition
u=max(x-Kx,zeros(N,1)); u0=u;

%% RBF
phi='gs';

%% Weights
W=sparse(N,N);

for ii=2:m
    se=x(ii);
    indc=1:n;
    sc=x(indc);
    
    wc=rbfga_weights_BS('BS1',ep,sc,se,{r,sig});
    W(ii,indc)=wc;
end

for ii=(m+1):(N-m)
    se=x(ii);
    indc=ii-m:ii+m;
    sc=x(indc);
    
    wc=rbfga_weights_BS('BS1',ep,sc,se,{r,sig});    
    W(ii,indc)=wc;
end

for ii=N-m+1:N-1
    se=x(ii);
    indc=N-n+1:N;
    sc=x(indc);
    
    wc=rbfga_weights_BS('BS1',ep,sc,se,{r,sig});
    W(ii,indc)=wc;
end

%% Integration
I=speye(N);

%BDF-1
A=I-W*dt;
[L,U]=lu(A);

u1=u;

b=u1;
b(end)=x(end)-Kx*exp(-r*dt);

u=U\(L\b);
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
indreg=[];
for ii=1:length(x)
    %         if (xfd(ii)-1)^2/((0.95*K)^2)+(yfd(ii)-1)^2/((0.95*K)^2)<=1
    if x(ii)>=1/3*Kx && x(ii)<=5/3*Kx
        indreg=[indreg ii];
    end
end

x=x(indreg);
u=u(indreg);

% K=100;
x=K*x;
u=K*u; u0=K*u0;

ua=rsol(sig, r, K, T, x);

% figure()
% plot(x,u0,'k--',x,u,'b-',x,ua,'r-')

err=(u-ua);

% figure()
% plot(x,abs(u-ua));


end

