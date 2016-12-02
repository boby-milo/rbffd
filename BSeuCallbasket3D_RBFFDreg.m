function [u,err,tim,x,dx,N,W] = BSeuCallbasket3D_RBFFDreg(Nx,n,ep,M)
%% 3D EU Call RBF-FD with BDF2

% Nx=30;
% ep=0.05;
tic
%% Parameters
% n=25;
% ep=2;
phi='gs';

% Nx=41;

%% Model
r=0.03;
sig1=0.15;
sig2=0.15;
sig3=0.15;
rho=[1      0.5     0.5;
    0.5     1       0.5;
    0.5     0.5     1];

T=1;
K=1;

%% Grid
display('Grid start')
Kx=1;

% Space
x=transpose(linspace(0,12*Kx,Nx));
dx=x(2)-x(1);
y=x;
z=x;

[X,Y,Z]=meshgrid(x,y,z);

xvec=transpose(X(:));
yvec=transpose(Y(:));
zvec=transpose(Z(:));


% ball=sqrt(xvec.^2+yvec.^2+zvec.^2);
ball=xvec+yvec+zvec;
ind=find(ball<=12+0.5*dx);

xvec=xvec(ind);
yvec=yvec(ind);
zvec=zvec(ind);
ball=ball(ind);

ind=1:numel(xvec);
indcf=1;

indff=find(ball>=max(ball)-0.99*dx); indff=unique(indff);
% indff=[find(xvec==max(xvec)),find(yvec==max(yvec)),find(zvec==max(zvec))]; indff=unique(indff);
indin=ind; indin([indff,indcf])=[];

N=numel(xvec);

tri=delaunay(xvec(indff),yvec(indff),zvec(indff));
figure(1)
plot3(xvec,yvec,zvec,'.')
hold on
axis equal
plot3(xvec(indff),yvec(indff),zvec(indff),'o')
plot3(xvec(indcf),yvec(indcf),zvec(indcf),'d')
grid on
plot3(xvec(indin),yvec(indin),zvec(indin),'o')

% pause()
% Time
dt=T/(M-1);
t=T:-dt:0;

display('Grid complete!')
%% Initial condition
u0=max((1/3)*(xvec+yvec+zvec)-Kx,zeros(1,length(xvec)));
u=u0';

%% RBF
s=[xvec' yvec' zvec'];
display('Distances')
indc=findKNearestNeighbors(s,s,n);

iind=repmat(indin,n,1); iind=iind(:); %n*N
jind=transpose(indc(indin,:)); jind=jind(:);%n*N
Wval=zeros(n,numel(indin));  %n*N

% Rc=xcdist(s,s,1);
% A=RBFmat(phi,ep,Rc,'0',1); display('A done')
%
%
% Ax=RBFmat(phi,ep,Rc,'1',1); display('Ax done')
% Ay=RBFmat(phi,ep,Rc,'1',2); display('Ay done')
% Az=RBFmat(phi,ep,Rc,'1',3); display('Az done')
%
% Axx=RBFmat(phi,ep,Rc,'2',1); display('Axx done')
% Ayy=RBFmat(phi,ep,Rc,'2',2); display('Ayy done')
% Azz=RBFmat(phi,ep,Rc,'2',3); display('Azz done')
%
% Axy=RBFmat(phi,ep,Rc,'m2',[1,2]); display('Axy done')
% Axz=RBFmat(phi,ep,Rc,'m2',[1,3]); display('Axz done')
% Ayz=RBFmat(phi,ep,Rc,'m2',[2,3]); display('Ayz done')

display('Stencils')
bb=0;
% internal points {
for ii=indin
    bb=bb+1;
    
    sc=[xvec(ii),yvec(ii),zvec(ii)]; xc=sc(:,1); yc=sc(:,2); zc=sc(:,3);
    se=s(indc(ii,:),:);
    
    Rc=xcdist(se,se,1);
    A=RBFmat(phi,ep,Rc,'0',1);
    
    Ax=RBFmat(phi,ep,Rc,'1',1); %display('Ax done')
    Ay=RBFmat(phi,ep,Rc,'1',2); %display('Ay done')
    Az=RBFmat(phi,ep,Rc,'1',3); %display('Az done')
    
    Axx=RBFmat(phi,ep,Rc,'2',1); %display('Axx done')
    Ayy=RBFmat(phi,ep,Rc,'2',2); %display('Ayy done')
    Azz=RBFmat(phi,ep,Rc,'2',3); %display('Azz done')
    
    Axy=RBFmat(phi,ep,Rc,'m2',[1,2]); %display('Axy done')
    Axz=RBFmat(phi,ep,Rc,'m2',[1,3]); %display('Axz done')
    Ayz=RBFmat(phi,ep,Rc,'m2',[2,3]); %display('Ayz done')
    
    l=transpose(-r*A(1,:)...
    +r*xc'.*Ax(1,:)+r*yc'.*Ay(1,:)+r*zc'.*Az(1,:)...
    +0.5*sig1.^2*xc'.^2.*Axx(1,:)...
    +0.5*sig2.^2*yc'.^2.*Ayy(1,:)...
    +0.5*sig3.^2*zc'.^2.*Azz(1,:)...
    +rho(1,2)*sig1*sig2*xc'.*yc'.*Axy(1,:)...
    +rho(1,3)*sig1*sig3*xc'.*zc'.*Axz(1,:)...
    +rho(2,3)*sig2*sig3*yc'.*zc'.*Ayz(1,:));
    
    wc=A\l;
%     wc=rbffd3(ii,xc,yc,zc,indc(ii,:),r,sig1,sig2,sig3,rho,A,Ax,Ay,Az,Axx,Ayy,Azz,Axy,Axz,Ayz);
    Wval(:,bb)=wc;
end
% } internal points
Wval=Wval(:);
W=sparse(iind,jind,Wval,N,N);
I=speye(size(W));
display('Wights complete!')

%% Integration
display('Integration start')
display('BDF1')
% BDF-1
u1=u;
A=I-dt*W;
% [L1, U1]=lu(A);

b=u1;
b(indff)=(1/3)*(xvec(indff)+yvec(indff)+zvec(indff))-Kx*exp(-r*dt);

% u=L1\b;
% u=U1\u;
u=A\b;
u=max(u,zeros(size(u)));

display('BDF2')
% BDF-2
A=I-(2/3)*dt*W;
rcm=symrcm(A);
A=A(rcm,rcm);
[L1, U1]=lu(A);
for ii=3:M
    u2=u1;
    u1=u;
    b=(4/3)*u1-(1/3)*u2;
    b(indff)=(1/3)*(xvec(indff)+yvec(indff)+zvec(indff))-Kx*exp(-r*(ii-1)*dt);
    
    u(rcm)=L1\b(rcm);
    u(rcm)=U1\u(rcm);
    
    u=max(u,zeros(size(u)));
end
tim=toc;
display('Integration complete!')

%% Visualization
U = griddata(xvec,yvec,zvec,u,X,Y,Z);

xslice=[0];
yslice=xslice;
zslice=xslice;

figure(2)
slice(X,Y,Z,U,xslice,yslice,zslice)
axis equal
% colormap jet
view(128,26)
err=0;
end
