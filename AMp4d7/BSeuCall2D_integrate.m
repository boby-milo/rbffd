function [u,xvec,yvec,Acond] = BSeuCall2D_integrate(u,I,T,M,W,xscale,yscale,xvec,yvec,indin,indff,Kx,r,method)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[k,beta0,beta1,beta2]=BDF2coeffs(T,M);

t = k(1);
u1 = u;

A = I - beta0*W;
Acond = condest(A);

switch method
    case 'backslash'
        rcm = symrcm(A);
        A = A(rcm,rcm);
        [L, U] = lu(A);

        b = u1;
        b(indff) = 0.5*(xvec(indff)+yvec(indff)) - Kx*exp(-r*t);
%         b(inddy0) = 2*u(inddy0-1) - u(inddy0-2);
%         b(inddy0) = xvec(inddy0);

        for ii = 1:M
%             [ii M]
            u(rcm) = L\b(rcm);
            u(rcm) = U\u(rcm);

            % prepare next time step
            jj = min(M,ii+1);
            t = t + k(jj);

            u2 = u1;
            u1 = u;

            b = beta1(jj)*u1 - beta2(jj)*u2;
            b(indff) = 0.5*(xvec(indff)+yvec(indff)) - Kx*exp(-r*t);
%             b(inddy0) = 2*u(inddy0-1) - u(inddy0-2);
%             b(inddy0) = xvec(inddy0);
        end

    case 'gmres'
        setup.type = 'nofill';
        [L, U] = ilu(A,setup);

        b = u1;
        b(indff) = 0.5*(xvec(indff)+yvec(indff)) - Kx*exp(-r*t);
%         b(inddy0) = 2*u(inddy0-1) - u(inddy0-2);

        for ii = 1:M
%             [ii M]
            [u,flag,~,~,~] = gmres(A,b,[],1e-6,200,L,U,u1);
            if flag >0
                pause();
            end
            u = max(u,0);
            % prepare next time step
            jj = min(M,ii+1);
            t = t + k(jj);

            u2 = u1;
            u1 = u;

            b = beta1(jj)*u1 - beta2(jj)*u2;
            b(indff) = 0.5*(xvec(indff)+yvec(indff)) - Kx*exp(-r*t);
%             b(inddy0) = 2*u(inddy0-1) - u(inddy0-2);
        end
end

xvec = xscale*xvec;
yvec = yscale*yvec;
u = xscale*u;

end
