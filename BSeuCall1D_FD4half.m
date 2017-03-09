function [u,err,tim,x,dx,n,N,W] = BSeuCall1D_FD4half(Nx,M,Kmul)
%K,T,r,sig,M
%% 1D European Call FDo4
% Copyright 2016, Slobodan Milovanovic
% 2016-07-07
funname = 'BSeuCall1D_FD4half';
datafilename = [funname,'___Nx',num2str(Nx),'_M',num2str(M),'_Kmul',num2str(Kmul),'.mat'];

names = dir(['./Data/',datafilename]);

if ~isempty(names)
    cd('./Data')
    load(datafilename);
    cd('..')
else
    tic
    %% Parameters
    K=1; %strike
    T=1; %maturation
    r=0.03; %interest
    sig=0.15; %volatility
    
    %% Grid
    % mul=4;
    x=transpose(linspace(0,Kmul*K,Nx));
    dx=x(2)-x(1);
    
    N = numel(x);
    
    dt=T/M;
    % t=T:-dt:0;
    
    %% Initial Conditions
    % u=max(x-K,zeros(N,1));
    
    u=rsol(sig,r,K,T/2,x);
    t=T/2;
    %% FD
    
    n = 5;
    
    aa = -(1/12)*(r*x)/dx +(1/24)*(sig^2*x.^2)/dx^2;
    bb = (2/3)*(r*x)/dx -(2/3)*(sig^2*x.^2)/dx^2;
    cc = (5/4)*(sig^2*x.^2)/dx^2 +r;
    dd = -(2/3)*(r*x)/dx -(2/3)*(sig^2*x.^2)/dx^2;
    ee = (1/12)*(r*x)/dx +(1/24)*(sig^2*x.^2)/dx^2;
    
    abcde = [[aa(3:end);0;0],[bb(2:end);0],cc,[0;dd(1:end-1)],[0;0;ee(1:end-2)]];
    d = [-2,-1,0,1,2];
    
    W = (spdiags(abcde,d,N,N));
    
    aal = (1/4)*(r*x(2))/dx -(11/12)*(sig^2*x(2).^2)/dx^2;
    bbl = (5/6)*(r*x(2))/dx +(5/3)*(sig^2*x(2).^2)/dx^2;
    ccl = -(3/2)*(r*x(2))/dx -(1/2)*(sig^2*x(2).^2)/dx^2 +r;
    ddl = (1/2)*(r*x(2))/dx -(1/3)*(sig^2*x(2).^2)/dx^2;
    eel = -(1/12)*(r*x(2))/dx +(1/12)*(sig^2*x(2).^2)/dx^2;
    
    aar = (1/12)*(r*x(end-1))/dx +(1/12)*(sig^2*x(end-1).^2)/dx^2;
    bbr = -(1/2)*(r*x(end-1))/dx -(1/3)*(sig^2*x(end-1).^2)/dx^2;
    ccr = (3/2)*(r*x(end-1))/dx -(1/2)*(sig^2*x(end-1).^2)/dx^2 +r;
    ddr = -(5/6)*(r*x(end-1))/dx +(5/3)*(sig^2*x(end-1).^2)/dx^2;
    eer = -(1/4)*(r*x(end-1))/dx -(11/12)*(sig^2*x(end-1).^2)/dx^2;
    
    W(1,:)=zeros(1,N); %W(1,1)=1;
    W(2,:)=zeros(1,N); W(2,1:5)=[aal,bbl,ccl,ddl,eel];
    W(end-1,:)=zeros(1,N); W(end-1,(end-4):end)=[aar,bbr,ccr,ddr,eer];
    W(end,:)=zeros(1,N); %W(end,end)=1;
    
    W=-W;
    
    %% Integration
    I=speye(N);
    
    %BDF-1
    t=t+dt;
    A=I-W*dt;
    u1=u;
    
    b=u1;
    ii=M/2+1;
    b(end)=x(end)-K*exp(-r*(ii-1)*dt);
    
    u=A\b;
    u=max(u,0);
    
    %BDF-2
    A=I-(2/3)*dt*W;
    rcm=symrcm(A);
    A=A(rcm,rcm);
    [L1, U1]=lu(A);
    for ii=((M/2)+2):M
        t=t+dt;
        u2=u1;
        u1=u;
        
        b=((4/3)*u1-(1/3)*u2);
        b(end)=x(end)-K*exp(-r*(ii-1)*dt);
        
        u(rcm)=L1\b(rcm);
        u(rcm)=U1\u(rcm);
        u=max(u,0);
    end
    tim=toc;
    
    %% Error
    ua=rsol(sig, r, K, T, x);
    
    % figure()
    % plot(x,u0,'k--',x,u,'b-',x,ua,'r-')
    
    err=u-ua;
    
    cd('./Data')
    save(datafilename, 'u', 'err', 'tim', 'x', 'dx', 'n', 'N', 'W', 'Nx', 'M', 'Kmul')
    cd('..')
end
end