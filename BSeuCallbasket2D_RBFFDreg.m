function [u,err,tim,x,dx,N,W] = BSeuCallbasket2D_RBFFDreg(Nx,n,ep,M)
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
Kx=1;
x=transpose(linspace(0,8*Kx,Nx));
% dx=x(2)-x(1)
y=x;

xvec=[]; yvec=[];
for ii=1:Nx
    xl=linspace(0,x(ii),ii);
    yl=linspace(x(ii),0,ii);

    xvec=[xvec,xl];
    yvec=[yvec,yl];
end

N=numel(xvec);

ind=1:N;
indcf=1;
indff=[length(xvec)-Nx+1:length(xvec)];
indin=ind; indin([indff,indcf])=[];

L=8;
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
% pause()

%% RBF
phi='gs';
s = [xvec' yvec'];
% Weights
indc=findKNearestNeighbors(s,s,n);

iind=repmat(indin,n,1); iind=iind(:); %n*N
jind=transpose(indc(indin,:)); jind=jind(:);%n*N
Wval=zeros(n,numel(indin));  %n*N

lc=zeros(n+1,1);

for ii=indin
    sc=s(ii,:); xc=sc(:,1); yc=sc(:,2);
    se=s(indc(ii,:),:);

    Rc=xcdist(se,se,1);
    A=RBFmat(phi,ep,Rc,'0',1);

    Ax=RBFmat(phi,ep,Rc,'1',1);
    Ay=RBFmat(phi,ep,Rc,'1',2);

    Axx=RBFmat(phi,ep,Rc,'2',1);
    Ayy=RBFmat(phi,ep,Rc,'2',2);
    Axy=RBFmat(phi,ep,Rc,'m2',1:2);
    
    o=ones(1,n);
    Ac=[A, transpose(o);
        o, 0];
    
    lc(1:n,1)=transpose(-r*A(1,:)...
        +r*xc'.*Ax(1,:)+r*yc'.*Ay(1,:)...
        +0.5*sig1^2*xc'.^2.*Axx(1,:)...
        +0.5*sig2^2*yc'.^2.*Ayy(1,:)...
        +rho*sig1*sig2*xc'.*yc'.*Axy(1,:));
    
    lc(n+1,1)=-r;

    wc=Ac\lc;
    Wval(:,ii-1)=wc(1:end-1);
end

Wval=Wval(:);
W=sparse(iind,jind,Wval,N,N);
I = speye(size(W));

%% Integration
% BDF-1
u1=u;
A=I-dt*W;

b=u1;
b(indff)=0.5*(xvec(indff)+yvec(indff))-Kx*exp(-r*dt);

u=A\b;
u=max(u,zeros(size(u)));

% BDF-2
A=I-(2/3)*dt*W;
rcm=symrcm(A);
A=A(rcm,rcm);
[L1, U1]=lu(A);
for ii=3:M
    u2=u1;
    u1=u;
    b=(4/3)*u1-(1/3)*u2;
    b(indff)=0.5*(xvec(indff)+yvec(indff))-Kx*exp(-r*(ii-1)*dt);

    u(rcm)=L1\b(rcm);
    u(rcm)=U1\u(rcm);

    u=max(u,zeros(size(u)));
end
tim=toc;


%% Error
indreg=[];
for ii=1:length(xvec)
    %         if (xfd(ii)-1)^2/((0.95*K)^2)+(yfd(ii)-1)^2/((0.95*K)^2)<=1
    if xvec(ii)>=1/3*Kx && xvec(ii)<=5/3*Kx && yvec(ii)>=1/3*Kx && yvec(ii)<=5/3*Kx
        indreg=[indreg ii];
    end
end

xvec=xvec(indreg);
yvec=yvec(indreg);
u=u(indreg);

x=[xvec' yvec'];


uinterp=griddata(xulti,yulti,uulti,xvec,yvec,'cubic');

err=uinterp'-u;

end
