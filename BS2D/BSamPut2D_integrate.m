function [u,xvec,yvec,Acond] = BSamPut2D_integrate(u,lambda,I,T,M,N,W,xscale,yscale,xvec,yvec,indin,indff,Kx,r,method)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Integration
dt=T/(M-1);
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
Acond = condest(A);
% rcm=symrcm(A);
% A=A(rcm,rcm);
% [L1, U1]=lu(A);

setup.droptol = 1e-3;
[LA, UA]=ilu(A,setup);

for ii=3:M
%     waitbar(ii/M)
    u2=u1;
    u1=u;

    b=(4/3)*u1 - (1/3)*u2 + (2/3)*dt*lambda;
    
%     util(rcm)=L1\b(rcm);
%     util(rcm)=U1\util(rcm);
    
    [util,~] = gmres(A,b,6,1e-6,200,LA,UA);
    
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

% [k,beta0,beta1,beta2]=BDF2coeffs(T,M);
% 
% t = k(1);
% u1 = u;
% 
% A = I - beta0*W;
% Acond = condest(A)
% 
% switch method
%     case 'backslash'
%         %         rcm = symrcm(A);
%         %         A = A(rcm,rcm);
%         %         [L, U] = lu(A);
%         %
%         %         b = u1;
%         %
%         %         for ii = 1:M
%         % %             [ii M]
%         %             u(rcm) = L\b(rcm);
%         %             u(rcm) = U\u(rcm);
%         %
%         %             util=A\(u1+dt*lambda);
%         %             lambdaold=lambda;
%         %             lambda=zeros(N,1);
%         %             u=util+dt*(lambda-lambdaold);
%         %
%         %
%         %
%         %             % prepare next time step
%         %             jj = min(M,ii+1);
%         %             t = t + k(jj);
%         %
%         %             u2 = u1;
%         %             u1 = u;
%         %
%         %             b = beta1(jj)*u1 - beta2(jj)*u2;
%         %             b(indff) = 0.5*(xvec(indff)+yvec(indff)) - Kx*exp(-r*t);
%         % %             b(inddy0) = 2*u(inddy0-1) - u(inddy0-2);
%         % %             b(inddy0) = xvec(inddy0);
%         %           end
%         
%     case 'gmres'
%         setup.type = 'nofill';
%         [L, U] = ilu(A,setup);
%         b = u1;
%         
%         for ii = 1:M
%             [util,flag,~,~,~] = gmres(A,b,[],1e-6,200,L,U,u1);
%             if flag >0
%                 pause();
%             end
%             lambdaold=lambda;
%             lambda=zeros(numel(u),1);
%             u=util+(2/3)*beta0*(lambda-lambdaold);
%             for jj=1:numel(u)
%                 if u(jj)-(Kx-0.5*(xvec(jj)+yvec(jj)))<0
%                     u(jj)=Kx-0.5*(xvec(jj)+yvec(jj));
%                     lambda(jj)=lambdaold(jj)+(3/(2*beta0))*(u(jj)-util(jj));
%                 end
%             end
%             
%             u = max(u,0);
%             % prepare next time step
%             jj = min(M,ii+1);
%             t = t + k(jj);
%             
%             u2 = u1;
%             u1 = u;
%             
%             b = beta1(jj)*u1 - beta2(jj)*u2;
%         end
% end

xvec = xscale*xvec;
yvec = yscale*yvec;
u = xscale*u;

end