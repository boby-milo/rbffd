function [u,err,tim,x,dx,n,N,W]=BSeuCall2Dbasket_FD(Nx,M,Kmul)
%% 2D EU Call FD with BDF2
% 05-Jul-2016
load('UrefEU.mat')

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
u0=max(0.5*(xvec'+yvec')-Kx,zeros(1,length(xvec)));
u=u0';

%% FD
n = 9;
aa=0.5*sig1^2*xvec.^2/dx^2;
aaa=r*xvec/(2*dx);
bb=0.5*sig2^2*yvec.^2/dy^2;
bbb=r*yvec/(2*dy);
cc=rho*sig1*sig2*xvec.*yvec/(4*dx*dy);

W=sparse(Nx^2);
B=[-cc, -aa+aaa, cc, -bb+bbb, 2*aa+2*bb+r, -bb-bbb, cc, -aa-aaa, -cc];

C=[[B(Nx+2:end,1); zeros(Nx+1,1)], [B(Nx+1:end,2); zeros(Nx,1)], [B(Nx:end,3); zeros(Nx-1,1)],...
    [B(2:end,4); zeros(1,1)], B(:,5), [zeros(1,1); B(1:end-1,6)],...
    [zeros(Nx-1,1); B(1:end-Nx+1,7)],[zeros(Nx,1); B(1:end-Nx,8)],[zeros(Nx+1,1); B(1:end-Nx-1,9)]];

diagind=[-Nx-1, -Nx, -Nx+1, -1, 0, 1, Nx-1, Nx, Nx+1];
W=spdiags(C,diagind,Nx^2,Nx^2);

%% Integration

% BDF-1
u1=u;

W(indff,:)=[]; W(:,indff)=[];
I=speye(size(W));
A=I+dt*W;
A(1,:)=zeros(1,(Nx-1)^2); A(1,1)=1;
% [LA UA]=lu(A);

%Far-field BC
u(indff)=0.5*(xvec(indff)+yvec(indff))-Kx*exp(-r*dt);

uin=u; uin(indff)=[];

bin=zeros(1,Nx^2);
bin(Nx-1)=-u(Nx) * (bb(Nx-1)+bbb(Nx-1));
bin(Nx^2-2*Nx+1)=-u(Nx^2-Nx+1) * (aa(Nx^2-2*Nx+1)+aaa(Nx^2-2*Nx+1));

indupin=indup(2:end-2)-1;
bin(indupin)=u(indupin+1-Nx).*cc(indupin) - u(indupin+1).*(bbb(indupin)+bb(indupin)) - u(indupin+1+Nx).*cc(indupin);

indriin=indri(2:end-2)-Nx;
bin(indriin)=-u(indriin+1+Nx).*cc(indriin) - u(indriin+Nx).*(aaa(indriin)+aa(indriin)) + u(indriin+Nx-1).*cc(indriin);

indffin=Nx^2-Nx-1;
bin(indffin)=u(indffin+1-Nx).*cc(indffin) - u(indffin+1).*(bbb(indffin)+bb(indffin)) - u(indffin+1+Nx).*cc(indffin)...
    - u(indffin+Nx).*(aaa(indffin)+aa(indffin)) + u(indffin+Nx-1).*cc(indffin);

b=transpose(u'-bin*dt);
b(indff)=[];

%Solve
% uin=UA\(LA\b);
uin=A\b;
U=reshape(uin,Nx-1,Nx-1);

U=[U, u(indri(1:end-1))];
U=[U; transpose(u(indup));];
u=reshape(U,Nx^2,1);
u=max(u,zeros(size(u)));

% BDF-2
A=I+(2/3)*dt*W;
[LA, UA]=lu(A);
for ii=3:M
    u2=u1; u2in=u2; u2in(indff)=[];
    u1=u;  u1in=u1; u1in(indff)=[];
    
    %Far-field BC
    u(indff)=0.5*(xvec(indff)+yvec(indff))-Kx*exp(-r*(ii-1)*dt);
    
    uin=u; uin(indff)=[];
    
    bin=zeros(1,Nx^2);
    bin(Nx-1)=-u(Nx) * (bb(Nx-1)+bbb(Nx-1));
    bin(Nx^2-2*Nx+1)=-u(Nx^2-Nx+1) * (aa(Nx^2-2*Nx+1)+aaa(Nx^2-2*Nx+1));
    
    indupin=indup(2:end-2)-1;
    bin(indupin)=u(indupin+1-Nx).*cc(indupin) - u(indupin+1).*(bbb(indupin)+bb(indupin)) - u(indupin+1+Nx).*cc(indupin);
    
    indriin=indri(2:end-2)-Nx;
    bin(indriin)=-u(indriin+1+Nx).*cc(indriin) - u(indriin+Nx).*(aaa(indriin)+aa(indriin)) + u(indriin+Nx-1).*cc(indriin);
    
    indffin=Nx^2-Nx-1;
    bin(indffin)=u(indffin+1-Nx).*cc(indffin) - u(indffin+1).*(bbb(indffin)+bb(indffin)) - u(indffin+1+Nx).*cc(indffin)...
        - u(indffin+Nx).*(aaa(indffin)+aa(indffin)) + u(indffin+Nx-1).*cc(indffin);
    
    b=transpose((4/3)*u1'-(1/3)*u2'-(2/3)*bin*dt);
    b(indff)=[];
    
    uin=UA\(LA\b);
    U=reshape(uin,Nx-1,Nx-1);
    
    U=[U, u(indri(1:end-1))]; %transpose(u(indff(1:Nx-1)))];
    U=[U; transpose(u(indup));];
    u=reshape(U,Nx^2,1);
    u=max(u,zeros(size(u)));
end
tim=toc;

%% Error
% indreg=[];
% for ii=1:length(xvec)
%     %         if (xfd(ii)-1)^2/((0.95*Kx)^2)+(yfd(ii)-1)^2/((0.95*Kx)^2)<=1
%     if xvec(ii)>=1/3*Kx && xvec(ii)<=5/3*Kx && yvec(ii)>=1/3*Kx && yvec(ii)<=5/3*Kx
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

uinterp=griddata(xulti,yulti,uulti,xvec,yvec,'cubic');

err=uinterp-u;
end