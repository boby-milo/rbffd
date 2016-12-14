%% 1D European Call RBF-FD
% Copyright 2015, Slobodan Milovanovic
% 2016-02-06
function [x,dx,err]=BSeuCall1D_FDnoobsmooth(N,M,Kmul)
%% Parameters
K=1;
Kx=1; %strike
T=1; %maturation
r=0.03; %interest
sig=0.15; %volatility

%% Grid
%N=8001;
% mul=25;
x=transpose(linspace(0,Kmul*K,N));
dx=x(2)-x(1);

indin=2:N-1;
%M=50;
%M=round((T/2+dt)/dt);
n=3; %stencil size
m=round((n-1)/2);

% M=100000;
dt=T/(M-1);
t=T:-dt:0;
%dt=(T/2)/(M-1);
%t=(T/2):-dt:0;

%% Initial condition
%u=rsol(sig, r, Kx, T/2, x);
% u=max(x-Kx,zeros(N,1));
% u=max(u,0);
% u0=u;


Kx = 1;
fu = @(x) max(x-Kx,0);

u = smooth4(x,fu);

% u=@(x) max(x-Kx,0);
% 
% %4th order smoothing
% ff1=@(x) 2*(1/12)*(pi/2)*(((x+2).^3).*sign(x+2) -((x-2).^3).*sign(2-x)); %fourier(cos(2*w)/w^4, w, -x)
% ff2=@(x) 2*(1/12)*(pi/2)*(((x+3).^3).*sign(x+3) -((x-3).^3).*sign(3-x)); %fourier(cos(3*w)/w^4, w, -x)
% ff3=@(x) 2*(1/12)*(pi/2)*(((x+1).^3).*sign(x+1) -((x-1).^3).*sign(1-x)); %fourier(cos(w)/w^4, w, -x)
% 
% f1=@(x) (4*ff1(x))/(2*pi);
% f2=@(x) (-ff2(x)/3)/(2*pi);
% f3=@(x) ((14*pi*x.^3.*sign(x))/9)/(2*pi);
% f4=@(x) (-13*ff3(x))/(2*pi);
% 
% M4=@(x) f1(x)+f2(x)+f3(x)+f4(x);
% for ii=1:N
%     
% %     util1(ii)=(1/dx) * integral(@(s)u(x(ii)-s), -dx/2,+dx/2);
% %     util2(ii)=(1/dx) * integral(@(s)(1-abs(s/dx)).*u(x(ii)-s), -dx,+dx);
%     util4(ii)=(1/(dx)) * integral(@(s)M4(s/dx).*u(x(ii)-s), -3*dx,+3*dx);
% end
% 
% u=util4'; u0=u;

% plot(x,u0)
% pause

%From ass 2 in Computational Finance
s=x;
ds=dx;
sigma=sig;
gamma=1;

alpha=0.5*sigma^2*s.^2/(ds^2);
beta=r*s/ds;
l2=1/12*alpha-1/12*beta;
l1=-4/3*alpha+2/3*beta;
d=r+5/2*alpha;
upper1=-4/3*alpha-2/3*beta;
upper2=1/12*alpha+1/12*beta;

%Forward/Backward of first inner points
cd_scheme=diag(l2(3:end),-2)+diag(l1(2:end),-1)+diag(d)+diag(upper1(1:end-1),1)+diag(upper2(1:end-2),2);
cd_scheme(1,:)=zeros(1,N);
cd_scheme(2,:)=zeros(1,N);
%cd_scheme(2,:)=[0, r-15/4*alpha(2)+25/12*beta(2), 77/6*alpha(3)-4*beta(3), -107/6*alpha(4)+3*beta(4), 13*alpha(5)-4/3*beta(5), -61/12*alpha(6)+1/4*beta(6), 5/6*alpha(7), zeros(1,N-7)];
%cd_scheme(end-1,:)=[zeros(1,N-7), (5/6*alpha(end-6)), (-61/12*alpha(end-5)+1/4*beta(end-5)), (13*alpha(end-4)-4/3*beta(end-4)), (-107/6*alpha(end-3)+3*beta(end-3)), (77/6*alpha(end-2)-4*beta(end-2)), r+(-15/4*alpha(end-1)+25/12*beta(end-1)), 0];
cd_scheme(end-1,:)=zeros(1,N);
cd_scheme(end,:)=zeros(1,N);
W=cd_scheme;

%% Integration
I=speye(N);

%BDF-1
A=I+W*dt;
%A(1,1)=1;
%A(end,end)=1;
u1=u;

b=u1;

b(end)=x(end)-Kx*exp(-r*(dt));
b(end-1)=x(end-1)-Kx*exp(-r*(dt));

%b(end)=x(end)-Kx*exp(-r*(dt+T/2));
%b(end-1)=x(end-1)-Kx*exp(-r*(dt+T/2));

u=A\b;
u=max(u,0);

%BDF-2
A=I+(2/3)*dt*W;
rcm=symrcm(A);
A=A(rcm,rcm);
[L1, U1]=lu(A);
for ii=3:M
%     if mod(ii,500)==0
%         ii
%     end
    u2=u1;
    u1=u;
    
    b=((4/3)*u1-(1/3)*u2);
    
    b(end)=x(end)-Kx*exp(-r*((ii-1)*dt));
    b(end-1)=x(end-1)-Kx*exp(-r*((ii-1)*dt));

 %   b(end)=x(end)-Kx*exp(-r*((ii-1)*dt+T/2));
 %   b(end-1)=x(end-1)-Kx*exp(-r*((ii-1)*dt+T/2));

    u(rcm)=L1\b(rcm);
    u(rcm)=U1\u(rcm);
    u=max(u,0);
end
%tim=toc;
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

ua=rsol(sig, r, K, T, x);

err=u-ua;