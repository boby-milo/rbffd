function [u,err,tim,x,dx,N,W,t] = BSeuCall1D_FDo4smooth(N,M,mul)
%K,T,r,sig,M
%% 1D European Call FDo4
% Copyright 2016, Slobodan Milovanovic
% 2016-07-07

tic
%% Parameters
K=1; %strike
T=1; %maturation
r=0.03; %interest
sig=0.15; %volatility

%% Grid
% mul=4;
x=transpose(linspace(0,mul*K,N));
dx=x(2)-x(1);

dt=T/(M-1);
% t=T:-dt:0;

%% Initial condition
K=1;
u=@(x) max(x-K,0);
t=0;

%4th order smoothing
ff1=@(x) 2*(1/12)*(pi/2)*(((x+2).^3).*sign(x+2) -((x-2).^3).*sign(2-x)); %fourier(cos(2*w)/w^4, w, -x)
ff2=@(x) 2*(1/12)*(pi/2)*(((x+3).^3).*sign(x+3) -((x-3).^3).*sign(3-x)); %fourier(cos(3*w)/w^4, w, -x)
ff3=@(x) 2*(1/12)*(pi/2)*(((x+1).^3).*sign(x+1) -((x-1).^3).*sign(1-x)); %fourier(cos(w)/w^4, w, -x)

f1=@(x) (4*ff1(x))/(2*pi);
f2=@(x) (-ff2(x)/3)/(2*pi);
f3=@(x) ((14*pi*x.^3.*sign(x))/9)/(2*pi);
f4=@(x) (-13*ff3(x))/(2*pi);

M4=@(x) f1(x)+f2(x)+f3(x)+f4(x);
for ii=1:N
    
%     util1(ii)=(1/dx) * integral(@(s)u(x(ii)-s), -dx/2,+dx/2);
%     util2(ii)=(1/dx) * integral(@(s)(1-abs(s/dx)).*u(x(ii)-s), -dx,+dx);
    util4(ii)=(1/(dx)) * integral(@(s)M4(s/dx).*u(x(ii)-s), -3*dx,+3*dx);
end

u=util4';

%% FD

% aa = 0.5*sig^2*x.^2/dx^2 -0.5*r*x/dx;
% bb = -sig^2*x.^2/dx^2 -r;
% cc = 0.5*sig^2*x.^2/dx^2 +0.5*r*x/dx;

aa = -(1/12)*(r*x)/dx +(1/24)*(sig^2*x.^2)/dx^2;
bb = (2/3)*(r*x)/dx -(2/3)*(sig^2*x.^2)/dx^2;
cc = (5/4)*(sig^2*x.^2)/dx^2 +r;
dd = -(2/3)*(r*x)/dx -(2/3)*(sig^2*x.^2)/dx^2;
ee = (1/12)*(r*x)/dx +(1/24)*(sig^2*x.^2)/dx^2;

abcde = [[aa(3:end);0;0],[bb(2:end);0],cc,[0;dd(1:end-1)],[0;0;ee(1:end-2)]];
d = [-2,-1,0,1,2];

W = (spdiags(abcde,d,N,N));

aal = (1/4)*(r*x(2))/dx -(11/12)*(sig^2*x(2).^2)/dx^2;
bbl = (5/6)*(r*x(2))/dx +(5/3)*(sig^2*x(2).^2)/dx^2;
ccl = -(3/2)*(r*x(2))/dx -(1/2)*(sig^2*x(2).^2)/dx^2 +r;
ddl = (1/2)*(r*x(2))/dx -(1/3)*(sig^2*x(2).^2)/dx^2;
eel = -(1/12)*(r*x(2))/dx +(1/12)*(sig^2*x(2).^2)/dx^2;

aar = (1/12)*(r*x(end-1))/dx +(1/12)*(sig^2*x(end-1).^2)/dx^2;
bbr = -(1/2)*(r*x(end-1))/dx -(1/3)*(sig^2*x(end-1).^2)/dx^2;
ccr = (3/2)*(r*x(end-1))/dx -(1/2)*(sig^2*x(end-1).^2)/dx^2 +r;
ddr = -(5/6)*(r*x(end-1))/dx +(5/3)*(sig^2*x(end-1).^2)/dx^2;
eer = -(1/4)*(r*x(end-1))/dx -(11/12)*(sig^2*x(end-1).^2)/dx^2;

W(1,:)=zeros(1,N); %W(1,1)=1;
W(2,:)=zeros(1,N); W(2,1:5)=[aal,bbl,ccl,ddl,eel];
W(end-1,:)=zeros(1,N); W(end-1,(end-4):end)=[aar,bbr,ccr,ddr,eer];
W(end,:)=zeros(1,N); %W(end,end)=1;

W=-W;

%% Integration
I=speye(N);

%BDF-1
A=I-W*dt;
t=t+dt;
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
    t=t+dt;
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

% figure()
% plot(x,u0,'k--',x,u,'b-',x,ua,'r-')

err=(u-ua);