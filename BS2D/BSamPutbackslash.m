function [u] = BSamPutbackslash(xvec,yvec,N,indin,indcf,indff,W,dt,u,r,Kx,T,M)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

lambda=zeros(N,1);
I = speye(size(W));
% BDF-1
u1 = u;
Abig = I-dt*W;
A = Abig(indin,indin);

u1=u1(indin);
lambda = lambda(indin);
util=A\(u1+dt*lambda);
lambdaold=lambda;
lambda=zeros(numel(indin),1);
u=util+dt*(lambda-lambdaold);

for ii=1:numel(indin)
    if u(ii)-(Kx-0.5*(xvec(indin(ii))+yvec(indin(ii))))<0
        u(ii)=Kx-0.5*(xvec(indin(ii))+yvec(indin(ii)));
        lambda(ii)=lambdaold(ii)+(u(ii)-util(ii))/dt;
    end
end
u=max(u,zeros(size(u)));

% BDF-2
Abig=I-(2/3)*dt*W;
A=Abig(indin,indin);
rcm=symrcm(A);
A=A(rcm,rcm);
[L1, U1]=lu(A);

for ii=3:M
    u2=u1;
    u1=u;

    b=(4/3)*u1 - (1/3)*u2 + (2/3)*dt*lambda;
    
    util(rcm)=L1\b(rcm);
    util(rcm)=U1\util(rcm);
    lambdaold=lambda;
    lambda=zeros(numel(indin),1);
    
    u=util+(2/3)*dt*(lambda-lambdaold);
    
    for jj=1:numel(indin)
        if u(jj)-(Kx-0.5*(xvec(indin(jj))+yvec(indin(jj))))<0
            u(jj)=Kx-0.5*(xvec(indin(jj))+yvec(indin(jj)));
            lambda(jj)=lambdaold(jj)+(3/(2*dt))*(u(jj)-util(jj));
        end
    end
    
    u=max(u,zeros(size(u)));
end
ufin(indcf) = Kx;
ufin(indin) = u;
ufin(indff) = zeros(size(indff));
u = ufin';
end