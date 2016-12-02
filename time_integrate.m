function [u,m] = time_integrate(c,p,m,s)
%Integrates the equation in time.
%Inputs:
%   m.time.type: 'BDF';
%   m.time.order: 1, 2;

%Outputs:
%   u: solution at t;

dt = c.T/m.time.M;
m.time.dt = dt;

u = zeros(s.Ntot,1);

switch m.time.type
    case 'BDF'
        I = speye(s.N);
        A = I - s.W*dt;
        u1 = s.u0;
        
        switch m.time.order
            case 1
                rcm = symrcm(A);
                A = A(rcm,rcm);
                [L,U] = lu(A);
                
                for ii = 2 : m.time.M
                    
                    b = u1;
                    b(end) = s.x(end) - c.K*exp(-p.r*(ii-1)*dt);
                    
                    u(rcm) = L\b(rcm);
                    u(rcm) = U\u(rcm);
                    
                    u1 = u;
                    %u = max(u,0);
                end
                return;
                
            case 2
                b = u1;
                b(end) = s.x(end) - c.K*exp(-p.r*dt);
                
                u = A\b;
                %u = max(u,0);
                
                %BDF-2
                A = I - (2/3)*dt*s.W;
                rcm = symrcm(A);
                A = A(rcm,rcm);
                [L,U] = lu(A);
                for ii = 3 : m.time.M
                    u2 = u1;
                    u1 = u;
                    
                    b = (4/3)*u1 - (1/3)*u2;
                    b(end) = s.x(end) - c.K*exp(-p.r*(ii-1)*dt);
                    
                    u(rcm) = L\b(rcm);
                    u(rcm) = U\u(rcm);
                    %u=max(u,0);
                end
                return;
        end
end

end

