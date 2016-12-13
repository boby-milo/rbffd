function [u,method] = time_integrate(contract,parameter,method,grid)
%Integrates the equation in time.
%Inputs:
%   method.time.type: 'BDF';
%   method.time.order: 1, 2;

%Outputs:
%   u: solution at t;

dt = contract.T/method.time.M;
method.time.dt = dt;

u = zeros(grid.Ntot,1);

switch method.time.type
    case 'BDF'
        I = speye(grid.N);
        A = I - grid.W*dt;
        u1 = grid.u0;
        
        switch method.time.order
            case 1
                rcm = symrcm(A);
                A = A(rcm,rcm);
                [L,U] = lu(A);
                
                for ii = 2 : method.time.M
                    
                    b = u1;
                    b(end) = grid.x(end) - contract.K*exp(-parameter.r*(ii-1)*dt);
                    
                    u(rcm) = L\b(rcm);
                    u(rcm) = U\u(rcm);
                    
                    u1 = u;
                    %u = max(u,0);
                end
                return;
                
            case 2
                b = u1;
                b(end) = grid.x(end) - contract.K*exp(-parameter.r*dt);
                
                u = A\b;
                %u = max(u,0);
                
                %BDF-2
                A = I - (2/3)*dt*grid.W;
                rcm = symrcm(A);
                A = A(rcm,rcm);
                [L,U] = lu(A);
                for ii = 3 : method.time.M
                    u2 = u1;
                    u1 = u;
                    
                    b = (4/3)*u1 - (1/3)*u2;
                    b(end) = grid.x(end) - contract.K*exp(-parameter.r*(ii-1)*dt);
                    
                    u(rcm) = L\b(rcm);
                    u(rcm) = U\u(rcm);
                    %u=max(u,0);
                end
                return;
        end
end

end

