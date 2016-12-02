function [W] = differentiation_matrix(m,p,s)
%Discretizes spatial differential operator.
%Inputs:
%   m.space.type: 'FD', 'RBF-FD';
%   m.space.order: 2, 4;
%Outputs:
%   W: differentiation matrix;

W = [];

switch m.space.type
    
    case 'FD'
        if s.type ~= 'uniform'
            disp('FD is implemented only for uniform grids!')
            return;
        else
            switch s.dim
                case 1
                    
                    switch m.space.order
                        case 2
                            aa = 0.5*p.sig^2*s.x.^2/s.h^2 -0.5*p.r*s.x/s.h;
                            bb = -p.sig^2*s.x.^2/s.h^2 -p.r;
                            cc = 0.5*p.sig^2*s.x.^2/s.h^2 +0.5*p.r*s.x/s.h;
                            
                            W = gallery('tridiag',aa(2:end),bb,cc(1:end-1));
                            W(1,:)=zeros(1,s.Ntot);
                            W(end,:)=zeros(1,s.Ntot);
                            return;
                            
                        case 4
                            W = [];
                            return;
                    end
                    
                case 2
                    switch m.space.order
                        case 2
                            aa = 0.5*p.sig(1)^2*s.x(:,1).^2/s.h(1)^2;
                            aaa = p.r*s.x(:,1)/(2*s.h(1));
                            bb = 0.5*p.sig(2)^2*s.x(:,2).^2/s.h(2)^2;
                            bbb = p.r*s.x(:,2)/(2*s.h(2));
                            cc = p.rho*p.sig(1)*p.sig(2)*s.x(:,1).*s.x(:,2)/(4*s.h(1)*s.h(2));
                            
                            B = [ -cc, -aa + aaa, cc, -bb+bbb, 2*aa + 2*bb + p.r, -bb -bbb, cc, -aa -aaa, -cc];
                            
                            C = [[B(s.N+2:end,1); zeros(s.N+1,1)], [B(s.N+1:end,2); zeros(s.N,1)], [B(s.N:end,3); zeros(s.N-1,1)],...
                                [B(2:end,4); zeros(1,1)], B(:,5), [zeros(1,1); B(1:end-1,6)],...
                                [zeros(s.N-1,1); B(1:end-s.N+1,7)],[zeros(s.N,1); B(1:end-s.N,8)],[zeros(s.N+1,1); B(1:end-s.N-1,9)]];
                            
                            diagind = [-s.N-1, -s.N, -s.N+1, -1, 0, 1, s.N-1, s.N, s.N+1];
                            W = spdiags(C,diagind,s.Ntot,s.Ntot);
                            return;
                            
                        case 4
                            W = [];
                            return;
                    end
                    
                case 'RBF-FD'
                    W = [];
                    return;
            end
        end
end

end
