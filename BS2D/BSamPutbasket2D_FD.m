function [u,err,tim,x,dx,n,N,W]=BSamPutbasket2D_FD(Nx)
%% 2D AM Put FD with BDF2

load('UrefAM.mat')
tic
%% Model
r=0.03;
T=1;
K=100;
sig1=0.15;
sig2=0.15;
rho=0.5;

%% Grid
% Nx=1000;
Kmul = 8;
Kx = 1/Kmul;

x=linspace(0,1,Nx); dx=x(2)-x(1);
y=linspace(0,1,Nx); dy=y(2)-y(1);

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
u0=max(Kx-0.5*(xvec'+yvec'),zeros(1,length(xvec)));
u=u0';

%% FD
n = 9;
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
[u] = BSamPutbackslash(xvec,yvec,N,indin,indcf,indff,W,dt,u,r,Kx,T,M);
tim=toc;

%% Error
% indreg=[];
% for ii=1:length(xvec)
%     %         if (xfd(ii)-1)^2/((0.95*K)^2)+(yfd(ii)-1)^2/((0.95*K)^2)<=1
%     if xvec(ii)>=1/3*K && xvec(ii)<=5/3*K && yvec(ii)>=1/3*K && yvec(ii)<=5/3*K
%         indreg=[indreg ii];
%     end
% end
% 
% xvec=xvec(indreg);
% yvec=yvec(indreg);
% u=u(indreg);

xvec = K*Kmul*xvec;
yvec = K*Kmul*yvec;
x = [xvec yvec];
dx = K*Kmul*dx;
u = K*Kmul*u;

% uinterp=griddata(K*xulti,K*yulti,K*uulti,xvec,yvec,'cubic');

% err=uinterp-u;
err=[];

end