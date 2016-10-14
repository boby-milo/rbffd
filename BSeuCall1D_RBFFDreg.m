function [u,err,tim,x,dx,N,W] = BSeuCall1D_RBFFDreg(N,n,ep,M)
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
x=transpose(linspace(0,4,N));
dx=x(2)-x(1);

indin=2:N-1;

m=round((n-1)/2);

dt=T/(M-1);
% t=T:-dt:0;

%% Initial condition
u=max(x-Kx,zeros(N,1)); %u0=u;

%% RBF
phi='gs';

% Rc=xcdist(x,x,1);
% A=RBFmat(phi,ep,Rc,'0',1);
% Ax=RBFmat(phi,ep,Rc,'1',1);
% Axx=RBFmat(phi,ep,Rc,'2',1);

Rc=zeros(N,2*n-1,2);

index=1:n-1;
Rc(index,1:n+index-1,:)=xcdist(x(index),x(1:n+index-1),1);

index=n:N-n+1;
Rc(index,:,:)=xcdist(x(index),x(index-n+1:index+n-1),1);

index=N-n+2:N;
Rc(index,index-N+n:2*n-1,:)=xcdist(x(index),x(index-n+1:N),1);

A=RBFmat(phi,ep,Rc,'0',1);
Ax=RBFmat(phi,ep,Rc,'1',1);
Axx=RBFmat(phi,ep,Rc,'2',1);

A = repmat(A(n,:),[N 1]);
A = spdiags(A,-n+1:n-1,N,N);
Ax = repmat(Ax(n,:),[N 1]);
Ax = spdiags(Ax,-n+1:n-1,N,N);
Axx = repmat(Axx(n,:),[N 1]);
Axx = spdiags(Axx,-n+1:n-1,N,N);

%% Weights
iind=repmat(indin,n,1); iind=iind(:);
jind=zeros(n,N-2);
Wval=zeros(n,N-2);

lc=zeros(n+1,1);

for ii=2:m
    xc=x(ii);
    indc=1:n;
    
    o=ones(1,n);
    Ac=[A(indc,indc), transpose(o);
        o, 0];
    lc(1:n,1)=transpose(-r*A(ii,indc)+r*xc.*Ax(ii,indc)+0.5*xc.^2.*sig^2.*Axx(ii,indc));
    lc(n+1,1)=-r;
    wc=Ac\lc;
    Wval(:,ii-1)=wc(1:end-1);
    
    %     Ac=A(indc,indc);
    %     lc=transpose(-r*A(ii,indc)+r*xc.*Ax(ii,indc)+0.5*xc.^2.*sig^2.*Axx(ii,indc));
    %     wc=Ac\lc;
    %     Wval(:,ii-1)=wc;
    
    jind(:,ii-1)=indc';
end

for ii=(m+1):(N-m)
    xc=x(ii);
    indc=ii-m:ii+m;
    
    o=ones(1,n);
    Ac=[A(indc,indc), transpose(o);
        o, 0];
    lc(1:n,1)=transpose(-r*A(ii,indc)+r*xc.*Ax(ii,indc)+0.5*xc.^2.*sig^2.*Axx(ii,indc));
    lc(n+1,1)=-r;
    wc=Ac\lc;
    Wval(:,ii-1)=wc(1:end-1);
    
    %     Ac=A(indc,indc);
    %     lc=transpose(-r*A(ii,indc)+r*xc.*Ax(ii,indc)+0.5*xc.^2.*sig^2.*Axx(ii,indc));
    %     wc=Ac\lc;
    %     Wval(:,ii-1)=wc;
    
jind(:,ii-1)=indc';
end

for ii=(N-m+1):(N-1)
    xc=x(ii);
    indc=N-n+1:N;
    
    o=ones(1,n);
    Ac=[A(indc,indc), transpose(o);
        o, 0];
    lc(1:n,1)=transpose(-r*A(ii,indc)+r*xc.*Ax(ii,indc)+0.5*xc.^2.*sig^2.*Axx(ii,indc));
    lc(n+1,1)=-r;
    wc=Ac\lc;
    Wval(:,ii-1)=wc(1:end-1);
    
    %     Ac=A(indc,indc);
    %     lc=transpose(-r*A(ii,indc)+r*xc.*Ax(ii,indc)+0.5*xc.^2.*sig^2.*Axx(ii,indc));
    %     wc=Ac\lc;
    %     Wval(:,ii-1)=wc;
    
jind(:,ii-1)=indc';
end
jind=jind(:);
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