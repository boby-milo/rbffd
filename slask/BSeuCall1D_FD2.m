function [u,err,tim,x,dx,N,W] = BSeuCall1D_FD2(N,M,mul)
%K,T,r,sig,M
%% 1D European Call FD
% Copyright 2016, Slobodan Milovanovic
% 2016-07-07

tic
%% Parameters
K=1; %strike
T=1; %maturation
r=0.03; %interest
sig=0.15; %volatility

%% Grid
% mul=25;
x=transpose(linspace(0,mul*K,N));
dx=x(2)-x(1);

dt=T/(M-1);
t=T:-dt:0;

%% Initial Conditions
u=max(x-K,zeros(N,1));

%% FD

aa = 0.5*sig^2*x.^2/dx^2 -0.5*r*x/dx;
bb = -sig^2*x.^2/dx^2 -r;
cc = 0.5*sig^2*x.^2/dx^2 +0.5*r*x/dx;

W = gallery('tridiag',aa(2:end),bb,cc(1:end-1));
W(1,:)=zeros(1,N); %W(1,1)=1;
W(end,:)=zeros(1,N); %W(end,end)=1;

%% Integration
I=speye(N);

%BDF-1
A=I-W*dt;

u1=u;

b=u1;
b(end)=x(end)-K*exp(-r*dt);

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
    b(end)=x(end)-K*exp(-r*(ii-1)*dt);
    
    u(rcm)=L1\b(rcm);
    u(rcm)=U1\u(rcm);
    u=max(u,0);
end
tim=toc;

%% Error
ua=rsol(sig, r, K, T, x);
err=(u-ua);
end