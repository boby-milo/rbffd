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
            switch m.space.order
                case 2
                    aa = 0.5*p.sig^2*s.x.^2/s.h^2 -0.5*p.r*s.x/s.h;
                    bb = -p.sig^2*s.x.^2/s.h^2 -p.r;
                    cc = 0.5*p.sig^2*s.x.^2/s.h^2 +0.5*p.r*s.x/s.h;
                    
                    W = gallery('tridiag',aa(2:end),bb,cc(1:end-1));
                    W(1,:)=zeros(1,s.N);
                    W(end,:)=zeros(1,s.N);
                    return;
                case 4
                    W = [];
                    return;
            end
        end
    case 'RBF-FD'
        W = [];
        return;
end
end

