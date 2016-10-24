function [u,err,tim,x,dx,N,W] = BSeuCallbasket2D_RBFFDadapEpsavg(Nx,n,M,fit)
%% 2D EU Call RBF-FD with BDF2
% 2016-02-04 sparse
load('UrefEU.mat')

tic
%% Parameters
phi='gs';

switch n
    case 9
        gamma=0.011; gamma=fit*gamma;
    case 13
        gamma=0.0145; gamma=fit*gamma;
    case 25
        gamma=0.0450; gamma=fit*gamma;
end

%% Model
r=0.03;
sig1=0.15;
sig2=0.15;
rho=0.5;

T=1;
K=1;

%% Grid
Kx=1;

% Nx=100;
i=1:Nx;
Ki=2*Kx;
S=4*Ki;

g=3; %tune this! 1,2,3,4,5

c=2*Ki/g;

dxi=(1/Nx)*(asinh((S-Ki)/c)-asinh(-Ki/c));
xi=asinh(-Ki/c)+i*dxi;
x=[0, Ki+c*sinh(xi)];
y=zeros(numel(x),1);
Kind=-x(x<=Ki)+Ki;

xvec=[]; yvec=[];
for ii=1:numel(x)
    xl=linspace(0,x(ii),ii);
    yl=linspace(x(ii),0,ii);
    
    xvec=[xvec,xl];
    yvec=[yvec,yl];
end

N=numel(xvec);
ind=1:N;

indcf=1;
indff=(N-numel(x)+1):N;
indin=ind; indin([indff,indcf])=[];

L=8;
Nlinsq=sqrt(ceil(N*2-sqrt(N*2)));
dx=L/(Nlinsq-1);

% M=100000;
dt=T/(M-1);
t=T:-dt:0;

%% Initial condition
u0=max(0.5*(xvec+yvec)-Kx,zeros(1,length(xvec)));
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
s=[xvec' yvec'];

% Weights
indc=findKNearestNeighbors(s,s,n);

iind=repmat(indin,n,1); iind=iind(:); %n*N
jind=transpose(indc(indin,:)); jind=jind(:);%n*N
Wval=zeros(n,numel(indin));  %n*N

% internal points {
bb=0;
for ii=indin
    bb=bb+1;
    %     ii
    %     showsten(1,Nx,xvec,yvec,indc); pause()
    sc=[xvec(ii),yvec(ii)]; xc=sc(:,1); yc=sc(:,2);
    se=s(indc(ii,:),:);
    
    Rc=xcdist(se,se,1);
    
    H=Rc(:,:,1);
    havg=mean(H(H>0));
    ep=gamma/havg;
    
    A=RBFmat(phi,ep,Rc,'0',1);
    
    Ax=RBFmat(phi,ep,Rc,'1',1);
    Ay=RBFmat(phi,ep,Rc,'1',2);
    
    Axx=RBFmat(phi,ep,Rc,'2',1);
    Ayy=RBFmat(phi,ep,Rc,'2',2);
    Axy=RBFmat(phi,ep,Rc,'m2',1:2);
    
    o = ones(1,n);
    Ac = [A, transpose(o);
        o, 0];
    
    lc(1:n,1) = transpose(-r*A(1,:)...
        +r*xc'.*Ax(1,:)+r*yc'.*Ay(1,:)...
        +0.5*sig1^2*xc'.^2.*Axx(1,:)...
        +0.5*sig2^2*yc'.^2.*Ayy(1,:)...
        +rho*sig1*sig2*xc'.*yc'.*Axy(1,:));
    
    lc(n+1,1) = -r;
    
    wc = Ac\lc;
    Wval(:,ii-1) = wc(1:end-1);
end
% } internal points
Wval=Wval(:);
W=sparse(iind,jind,Wval,N,N);

%         display('Weights completed');
I=speye(size(W));

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
% figure(3) %solution
% tri = delaunay(xvec,yvec);
% trisurf(tri, xvec', yvec', u);
% axis vis3d
% axis tight
% xlabel('S_1');
% ylabel('S_2');
% zlabel('V(S_1,S_2)');
% drawnow;

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

x=[xvec' yvec'];
u=u(indreg);

uinterp=griddata(xulti,yulti,uulti,xvec,yvec,'cubic');

err=uinterp'-u;

% errorreg=max(abs(uinterp'-u));
% errornormreg=rms(uinterp'-u);

% display([ep,errorreg]);
% figure(4) %error plot
% tri = delaunay(xvec',yvec');
% trisurf(tri, xvec', yvec', abs(uinterp'-u));
% shading interp
% colorbar
% view(2)
% axis vis3d
% % axis equal
% axis tight
% xlabel('S_1');
% ylabel('S_2');
% zlabel('\Delta V(S_1,S_2)');
end
