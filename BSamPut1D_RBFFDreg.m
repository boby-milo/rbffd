function [u,err,tim,x,dx,N,W] = BSamPut1D_RBFFDreg(N,n,ep,M,parallel)
%% 1D American Put RBF-FD
% 2016-02-06

tic
%% Parameters
K=100;
Kx=1; %strike
T=1; %maturation
r=0.03; %interest
sig=0.15; %volatility

%% Grid
x=transpose(linspace(0,4,N));
dx=x(2)-x(1);

indin=2:N-1;

dt=T/(M-1);
% t=T:-dt:0;

%% Initial condition
u=max(Kx-x,zeros(N,1)); %u0=u;
lambda=zeros(N,1);

%% RBF
phi ='gs';
% parallel = 0;
W = BSweights1Drbffd(r,sig,x,N,n,indin,phi,ep,parallel);

%% Integration with Operator Splitting
I=speye(N);

%BDF-1
A=I-W*dt;

u1=u;
util=A\(u1+dt*lambda);

lambdaold=lambda;
u=util+dt*(lambda-lambdaold);

for ii=1:N
    if u(ii)-(Kx-x(ii))<0
        u(ii)=Kx-x(ii);
        lambda(ii)=lambdaold(ii)+(u(ii)-util(ii))/dt;
    end
end

u=max(u,0);

%BDF-2
A=I-(2/3)*dt*W;
rcm=symrcm(A);
A=A(rcm,rcm);
[L1, U1]=lu(A);
for ii=3:M
    u2=u1;
    u1=u;
    
    b=(4/3)*u1 - (1/3)*u2 + (2/3)*dt*lambda;
    
    util(rcm)=L1\b(rcm);
    util(rcm)=U1\util(rcm);
    lambdaold=lambda;
    lambda=zeros(N,1);
    
    u=util+(2/3)*dt*(lambda-lambdaold);
    
    
    for jj=1:N
        if u(jj)-(Kx-x(jj))<0
            u(jj)=Kx-x(jj);
            lambda(jj)=lambdaold(jj)+(3/(2*dt))*(u(jj)-util(jj));
        end
    end
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
u=K*u; %u0=K*u0;

% ua=rsol(sig, r, K, T, x); %% NOT CORRECT!

% figure()
% plot(x,u0,'k--',x,u,'b-',x,ua,'r-')

% err=(u-ua);
err=[];

% figure()
% plot(x,abs(u-ua));

end
