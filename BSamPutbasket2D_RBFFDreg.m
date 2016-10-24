function [u,err,tim,x,dx,N,W] = BSamPutbasket2D_RBFFDreg(Nx,n,ep,M)
%% 2D American Put RBF-FD with BDF2
% 2016-02-04 sparse
load('UrefAM.mat')

tic
%% Parameters
phi='gs';

%% Model
r=0.03;
sig1=0.15;
sig2=0.15;
rho=0.5;

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
u0=max(Kx-0.5*(xvec+yvec),zeros(1,length(xvec)));
u=u0';

lambda=zeros(N,1);

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
% 
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
s = [xvec' yvec'];
W = BSweights2Drbffd(r,sig1,sig2,rho,s,N,n,indin,phi,ep);
I = speye(size(W));

%% Integration
% BDF-1
u1=u;
A=I-dt*W;

util=A\(u1+dt*lambda);

lambdaold=lambda;
lambda=zeros(N,1);
u=util+dt*(lambda-lambdaold);

for ii=1:N
    if u(ii)-(Kx-0.5*(xvec(ii)+yvec(ii)))<0
        u(ii)=Kx-0.5*(xvec(ii)+yvec(ii));
        lambda(ii)=lambdaold(ii)+(u(ii)-util(ii))/dt;
    end
end

u=max(u,zeros(size(u)));

% BDF-2
A=I-(2/3)*dt*W;
rcm=symrcm(A);
A=A(rcm,rcm);
[L1, U1]=lu(A);
for ii=3:M
%     waitbar(ii/M)
    u2=u1;
    u1=u;

    b=(4/3)*u1 - (1/3)*u2 + (2/3)*dt*lambda;
    
    util(rcm)=L1\b(rcm);
    util(rcm)=U1\util(rcm);
    lambdaold=lambda;
    lambda=zeros(N,1);
    
    u=util+(2/3)*dt*(lambda-lambdaold);
    
    for jj=1:N
        if u(jj)-(Kx-0.5*(xvec(jj)+yvec(jj)))<0
            u(jj)=Kx-0.5*(xvec(jj)+yvec(jj));
            lambda(jj)=lambdaold(jj)+(3/(2*dt))*(u(jj)-util(jj));
        end
    end
    
    u=max(u,zeros(size(u)));
end
tim=toc;

% figure(3) %solution
% tri = delaunay(xvec,yvec);
% trisurf(tri, xvec', yvec', u);
% axis vis3d
% axis tight
% xlabel('s_1');
% ylabel('s_2');
% zlabel('u(s_1,s_2)');
% drawnow;
%% Error
% indreg=[];
% for ii=1:length(xvec)
%     %         if (xfd(ii)-1)^2/((0.95*K)^2)+(yfd(ii)-1)^2/((0.95*K)^2)<=1
%     if xvec(ii)>=1/3*Kx && xvec(ii)<=5/3*Kx && yvec(ii)>=1/3*Kx && yvec(ii)<=5/3*Kx
%         indreg=[indreg ii];
%     end
% end
% 
% xvec=xvec(indreg);
% yvec=yvec(indreg);
% u=u(indreg);

x=[xvec' yvec'];

uinterp=griddata(xulti,yulti,uulti,xvec,yvec,'cubic');

err=uinterp'-u;

% errorreg=max(abs(uinterp'-u));
% errornormreg=rms(uinterp'-u);

% display([ep,errorreg]);
% figure() %error plot
% tri = delaunay(xvec',yvec');
% trisurf(tri, xvec', yvec', abs(uinterp'-u));
% shading interp
% colorbar
% view(2)
% axis vis3d
% % axis equal
% axis tight
% xlabel('s_1');
% ylabel('s_2');
% zlabel('\Delta V(s_1,s_2)');
end