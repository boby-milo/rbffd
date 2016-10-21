function [u,err,tim,x,dx,N,W] = BSeuCall1D_RBFFDadapConst(Nx,n,ep,M,g)
%K,T,r,sig,M
%% 1D European Call RBF-FD
% Copyright 2015, Slobodan Milovanovic 
% 2016-02-06

tic
%% Parameters
K=100;
Kx=1; %strike
T=1; %maturation
r=0.03; %interest
sig=0.15; %volatility

%% Grid

% Nx=100;
i=1:Nx;
Ki=Kx;
S=4*Ki;

% g=5; %tune this! 1,2,3,4,5
c=2*Ki/g;

dxi=(1/Nx)*(asinh((S-Ki)/c)-asinh(-Ki/c));
xi=asinh(-Ki/c)+i*dxi;
x=[0, Ki+c*sinh(xi)]';

N=numel(x);
ind=1:N;
indcf=1;
indff=N;
indin=ind; indin([indff,indcf])=[];

L=4;
dx=L/(N-1);

% n=3; %stencil size
m=round((n-1)/2);

% M=100000;
dt=T/(M-1);
% t=T:-dt:0;

% figure()
% plot(x, zeros(size(x)),'.');
% axis tight

%% Initial condition
u=max(x-Kx,zeros(N,1)); %u0=u;

%% RBF
phi='gs';

% Rc=xcdist(x,x,1); 
% A=RBFmat(phi,ep,Rc,'0',1);
% Ax=RBFmat(phi,ep,Rc,'1',1);
% Axx=RBFmat(phi,ep,Rc,'2',1);

%% Weights
iind=repmat(indin,n,1); iind=iind(:);
jind=zeros((N-2)*n,1);
Wval=zeros(n,N-2);

lc=zeros(n+1,1);

bb=0;
for ii=2:m
    bb=bb+1;
    xc=x(ii);
    indc=1:n;
    
    Rc=xcdist(x(indc),x(indc),1);
    A=RBFmat(phi,ep,Rc,'0',1);
    Ax=RBFmat(phi,ep,Rc,'1',1);
    Axx=RBFmat(phi,ep,Rc,'2',1);
    
    
    o=ones(1,n);
    Ac=[A, transpose(o);
        o, 0];
    
    lc(1:n,1)=transpose(-r*A(ii,:)+r*xc.*Ax(ii,:)+0.5*xc.^2.*sig^2.*Axx(ii,:));
    lc(n+1,1)=-r;
    
    wc=Ac\lc;
    
    Wval(:,ii-1)=wc(1:end-1);
    jind(bb:bb+n-1)=indc;
    bb=bb+n-1;
end

for ii=(m+1):(N-m)
    bb=bb+1;
    xc=x(ii);
    indc=ii-m:ii+m;
    
    Rc=xcdist(x(indc),x(indc),1);
    
    A=RBFmat(phi,ep,Rc,'0',1);
    Ax=RBFmat(phi,ep,Rc,'1',1);
    Axx=RBFmat(phi,ep,Rc,'2',1);
    
    o=ones(1,n);
    Ac=[A, transpose(o);
        o, 0];
    
    lc(1:n,1)=transpose(-r*A(m+1,:)+r*xc.*Ax(m+1,:)+0.5*xc.^2.*sig^2.*Axx(m+1,:));
    lc(n+1,1)=-r;
    
    wc=Ac\lc;
    Wval(:,ii-1)=wc(1:end-1);
    jind(bb:bb+n-1)=indc;
    bb=bb+n-1;
end

for ii=(N-m+1):(N-1)
    bb=bb+1;
    xc=x(ii);
    indc=N-n+1:N;
    
    Rc=xcdist(x(indc),x(indc),1);
    
    A=RBFmat(phi,ep,Rc,'0',1);
    Ax=RBFmat(phi,ep,Rc,'1',1);
    Axx=RBFmat(phi,ep,Rc,'2',1);
    
    o=ones(1,n);
    Ac=[A, transpose(o);
        o, 0];
    
    lc(1:n,1)=transpose(-r*A(ii-N+n,:)+r*xc.*Ax(ii-N+n,:)+0.5*xc.^2.*sig^2.*Axx(ii-N+n,:));
    lc(n+1,1)=-r;
    
    wc=Ac\lc;
    Wval(:,ii-1)=wc(1:end-1);
    jind(bb:bb+n-1)=indc;
    bb=bb+n-1;
end

Wval=Wval(:);
W=sparse(iind,jind,Wval,N,N);


% neighbours=findKNearestNeighbors(x,x(2:N-1),n);
% for ii=2:N-1
%     indc=neighbours(ii-1,:);
%     xc=x(ii);
%  
% %     figure()
% %     clf
% %     plot(x,zeros(N,1),'bo');
% %     hold on
% %     plot(x(indc),zeros(length(indc),1),'r*');
% %     pause(0.01)
%     
%     Ac=A(indc,indc);
%     lc=transpose(-r*A(ii,indc)+r*xc.*Ax(ii,indc)+0.5*xc.^2.*sig^2.*Axx(ii,indc));
%     wc=Ac\lc;
%     
%     W(ii,indc)=wc;
% end

%% Integration
I=speye(N);

%BDF-1
A=I-W*dt;
u1=u;
b=u1;
b(end)=x(end)-Kx*exp(-r*dt);

u=A\b;
u=max(u,0);

%BDF-2
A=I-(2/3)*dt*W;
rcm=symrcm(A);
A=A(rcm,rcm);
[L1, U1]=lu(A);
for ii=3:M
    u2=u1;
    u1=u;

    b=((4/3)*u1-(1/3)*u2);
    b(end)=x(end)-Kx*exp(-r*(ii-1)*dt);

    u(rcm)=L1\b(rcm);
    u(rcm)=U1\u(rcm);
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

ua=rsol(sig, r, K, T, x);

% figure()
% plot(x,u0,'k--',x,u,'b-',x,ua,'r-')

err=(u-ua);

% figure()
% plot(x,abs(u-ua));


end

