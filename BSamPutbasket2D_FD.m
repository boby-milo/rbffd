function [u,err,tim,x,dx,N,W]=BSamPutbasket2D_FD(Nx,M)
%% 2D AM Put FD with BDF2
% 28-Mar-2016
load('UrefAM.mat')

tic
%% Model
r=0.03;
T=1;
K=1;
sig1=0.15;
sig2=0.15;
rho=0.5;

%% Grid
% Nx=1000;

x=linspace(0,8*K,Nx); dx=x(2)-x(1);
y=linspace(0,8*K,Nx); dy=y(2)-y(1);

[X,Y]=meshgrid(x,y);
xvec=X(:);
yvec=Y(:);

ind=1:numel(xvec);

indle=[1:Nx];
inddo=[1:Nx:Nx^2-Nx+1];
indcf=1;

indup=[Nx:Nx:Nx^2];
indri=[Nx^2-Nx+1:1:Nx^2];
indff=[indup(1:end-1),indri];

indin=ind; indin([indff,indcf])=[];

N=numel(xvec);

t=linspace(T,0,M);
dt=t(1)-t(2);

%% Initial condition
ii=1;
u0=max(K-0.5*(xvec'+yvec'),zeros(1,length(xvec)));
u=u0';

%% FD
aa=0.5*sig1^2*xvec.^2/dx^2;
aaa=r*xvec/(2*dx);
bb=0.5*sig2^2*yvec.^2/dy^2;
bbb=r*yvec/(2*dy);
cc=rho*sig1*sig2*xvec.*yvec/(4*dx*dy);

W=sparse(Nx^2,Nx^2);
B=[-cc, -aa+aaa, cc, -bb+bbb, 2*aa+2*bb+r, -bb-bbb, cc, -aa-aaa, -cc];

C=-[[B(Nx+2:end,1); zeros(Nx+1,1)], [B(Nx+1:end,2); zeros(Nx,1)], [B(Nx:end,3); zeros(Nx-1,1)],...
    [B(2:end,4); zeros(1,1)], B(:,5), [zeros(1,1); B(1:end-1,6)],...
    [zeros(Nx-1,1); B(1:end-Nx+1,7)],[zeros(Nx,1); B(1:end-Nx,8)],[zeros(Nx+1,1); B(1:end-Nx-1,9)]];

diagind=[-Nx-1, -Nx, -Nx+1, -1, 0, 1, Nx-1, Nx, Nx+1];
W=spdiags(C,diagind,Nx^2,Nx^2);

%% Integration
I=speye(size(W));
% BDF-1
u1=u;

A=I-dt*W;
A(1,:)=0; A(1,1)=1;
A(indff,:)=0;
A(:,indff)=0;
idx = sub2ind(size(A), indff, indff);
A(idx)=1;

lambda=zeros(N,1);

util=A\(u1-dt*lambda);
lambdaold=lambda;
lambda=zeros(N,1);
u=util+dt*(lambdaold-lambda);

for ii=1:N
    if u(ii)-(K-0.5*(xvec(ii)+yvec(ii)))<0
        u(ii)=K-0.5*(xvec(ii)+yvec(ii));
        lambda(ii)=lambdaold(ii)-(u(ii)-util(ii))/dt;
    end
end

u=max(u,zeros(size(u)));

% BDF-2
A=I-(2/3)*dt*W;
A(1,:)=0; A(1,1)=1;
A(indff,:)=0;
A(:,indff)=0;
idx = sub2ind(size(A), indff, indff);
A(idx)=1;

rcm=symrcm(A);
A=A(rcm,rcm);
[L1, U1]=lu(A);

for ii=3:M
    u2=u1;
    u1=u;
    
    b=(4/3)*u1-(1/3)*u2-(2/3)*dt*lambda;
    util(rcm)=L1\b(rcm);
    util(rcm)=U1\util(rcm);
    lambdaold=lambda;
    lambda=zeros(N,1);
    
    u=util+(2/3)*dt*(lambdaold-lambda);
    
    for jj=1:numel(lambda)
        if u(jj)-(K-0.5*(xvec(jj)+yvec(jj)))<0
            u(jj)=K-0.5*(xvec(jj)+yvec(jj));
            lambda(jj)=lambdaold(jj)+(3/(2*dt))*(util(jj)-u(jj));
        end
    end
    
    u=max(u,zeros(size(u)));
end
tim=toc;

%% Error
indreg=[];
for ii=1:length(xvec)
    %         if (xfd(ii)-1)^2/((0.95*K)^2)+(yfd(ii)-1)^2/((0.95*K)^2)<=1
    if xvec(ii)>=1/3*K && xvec(ii)<=5/3*K && yvec(ii)>=1/3*K && yvec(ii)<=5/3*K
        indreg=[indreg ii];
    end
end

xvec=xvec(indreg);
yvec=yvec(indreg);
u=u(indreg);

x=[xvec yvec];

uinterp=griddata(K*xulti,K*yulti,K*uulti,xvec,yvec,'cubic');

err=uinterp-u;
end